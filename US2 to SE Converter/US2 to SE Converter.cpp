#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include "dirent.h"
#include <algorithm> // for sort()
// #include <limits> // for max double

#define STARCUTOFF 3e28 // 2.5e28 ~ 13 Jupiter masses is the Brown Dwarf cutoff used by SE
#define HYDROUPPER 4e20 // Vesta < mass < Ceres of ~4e20 is used as definite upper cutoff
#define HYDROLOWER 2e19 // mass < Mimas of ~2e19 is used as lower cutoff depending on neighborhood
#define K_KGKGAU 7.53445e53
// In solar-mass, earth-mass & AU units, the PI discriminant constant K = 807.
// But we use kg, kg, AU instead, so our constant = K*pow(solar-mass, 5.0/2.0)/earth-mass = 7.53445e53

class StateVect
{
public:
	double x, y, z;

	StateVect(double a = 0.0, double b = 0.0, double c = 0.0)
	{
		x = a;
		y = b;
		z = c;
	}

	double magnitude()
	{
		return (sqrt(x * x + y * y + z * z));
	}
};

class Color
{
public:
	double r, g, b, a;

	Color(double x = -1.0, double y = 0.0, double z = 0.0, double w = 0.0)
	{
		r = x;
		g = y;
		b = z;
		a = w;
	}

	std::string Print()
	{
		return ("(" + std::to_string(r) + ", " + std::to_string(g) +", " + std::to_string(b) + ", " + std::to_string(a) + ")");
	}
};

struct Quaternion
{
	double w, x, y, z;
};

struct Object
{
	std::string name, type, class_;
	Object* parent;
	Object* partner; // for binaries and objects of similar mass, otherwise it's == parent
	std::vector<Object*> child;
	double mass, radius;

	// orbit stuff
	StateVect position, velocity, angularVelocity;
	Quaternion orientation;
	double semimajor, period, inclination, eccentricity, argOfPeriapsis, longOfAscNode, meanAnomaly, hillSphereRadius;
	double rotationPeriod, obliquity;

	// general
	int temp;
	double silicateMass, ironMass, waterMass, hydrogenMass, surfacePressure, albedo, greenhouse, roughness;

	// used for stars
	bool isStar;
	double luminosity;

	// color stuff
	Color colorBeach, colorLowland, colorUpland; // for rocky planets
	Color colorLayer0, colorLayer1, colorLayer2, colorLayer3, colorLayer4, colorLayer5, colorLayer6, colorLayer7; // for gass giants
	Color atmosphere;

	// in hydrostatic equilibrium?
	bool isRound;
};
std::string::size_type sz;
std::vector<Object> object;
Object* root;

double G; // gravitation constant imported from Universe Sandbox
const double PI = 3.141592653589793;

void GetData(std::ifstream&);
void BuildHierarchy(std::list<Object>&, std::list<Object>::iterator&);
void Bond(std::list<Object>&, std::list<Object>::iterator&, Object&, Object&);
void Typifier(Object&);
void Classifier(Object&);
void PrintFile(std::ofstream&, Object&);

// math functions
void UpdateBinary(Object&, Object&, Object&);
void CalcOrbit(Object&);
void CalcMoreOrbit(Object&);
bool ClearanceCheck(Object&);
StateVect To_EulerAngles(Quaternion q);
StateVect RotateVector(Quaternion&, StateVect&);
double Distance(Object&, Object&);
StateVect CrossProduct(StateVect&, StateVect&);
double DotProduct(StateVect&, StateVect&);
StateVect Scale(double, StateVect&);
StateVect Subtract(StateVect&, StateVect&);
StateVect Add(StateVect&, StateVect&);
StateVect Normalize(StateVect&);
double To_Degree(double&);

int main()
{
	std::string systemName, starFileName, planetFileName, holder;

	DIR *dir = opendir("input/");
	struct dirent *entry;
	if (dir)
	{
		// find all files in input
		while ((entry = readdir(dir)) != NULL)
		{
			// open files with .json extension
			systemName = entry->d_name;
			if (systemName.find(".json", systemName.size() - 5) == std::string::npos)
				continue;
			std::ifstream inputFile("input/" + systemName);
			if (!inputFile)
			{
				std::cout << "\nThere was an error opening a simulation file!\n";
				continue;
			}

			// finds the "Name": part of the simulation file and uses it for star and planet .sc files
			while (inputFile >> holder && holder != "},");
			while (inputFile >> holder && holder != "},");
			std::getline(inputFile, holder);
			std::getline(inputFile, holder);
			holder.erase(0, 8);
			holder.erase(holder.size() - 2, 2);
			// use file name if the simulation is unnamed
			if (holder != "Empty Universe") systemName = holder;
			else systemName.erase(systemName.size() - 5, 5);
			planetFileName = starFileName = systemName;
			starFileName = "output/" + starFileName + " Star.sc";
			planetFileName = "output/" + planetFileName + " Planet.sc";

			// finds gravitational constant
			while (inputFile >> holder && !(holder.find("\"Gravity\":") + 1));
			holder.erase(0, 10);
			G = std::stod(holder, &sz);

			// inputs every object into the global vector
			GetData(inputFile);
			inputFile.close();

			std::list<Object> binary;
			Object systemO;
			systemO.name = systemName + " System";
			binary.push_back(systemO);
			std::list<Object>::iterator b = binary.begin(); // binary object counter

			// this is where the magic happens
			BuildHierarchy(binary, b);
			// this is where the magic happens #2
			CalcOrbit(*root);
			// more magic requires binaries to already have some orbit info
			CalcMoreOrbit(*root);

			for (int i = 0; i < object.size(); i++)
			{
				Typifier(object.at(i));
				Classifier(object.at(i));
			}

			std::ofstream starFile(starFileName.c_str());
			starFile << "StarBarycenter\t\t\"" << binary.begin()->name << "\"\n{}\n";
			starFile.close();

			root->parent = &*binary.begin();
			std::ofstream planetFile(planetFileName.c_str());

			planetFile << "\n\n\t/*\tREAD ME!\n\n"
				<< "\tIf you are opening this file because you noticed that some of the objects from your simulation were converted wrong,\n"
				<< "\tyou will be able to fix them here. Because of the complicated nature of this conversion, it is impossible for every value\n"
				<< "\tto be exact, so if you REALLY want all the numbers in this file to be 100% the same as from your simulation, you will\n"
				<< "\tneed to go through every object and compare every value to the one you see in Universe Sandbox. However, most of the\n"
				<< "\tvalues are close, so I will tell you what values are most likely to be inaccurate.\n\n"

				<< "\tFirst off, if there were any binary objects in your simulation (stars orbiting each other, double planet systems, etc.), their\n"
				<< "\torbits will be the most inaccurate. You will need to find these objects (use ctrl + F) and look at their 'Orbits' tag. The first thing\n"
				<< "\tto fix will be their eccentricity. You can find this value easily in US. ArgOfPericenter and MeanAnomaly are probably incorrect as well.\n"
				<< "\tOf course, because these are binary objects, make sure their orbital parameters are EXACTLY the same (period, eccentricity, inclination, AscendingNode,\n"
				<< "\tand MeanAnomaly) and make sure one of the objects has an ArgOfPericenter +/- 180 degrees from the other.\n\n"

				<< "\tFor all other objects, there are some special cases where values are calculated incorrectly. If the Ascending Node of an object in US\n"
				<< "\tis undefined, then you can ignore the value printed in this file, but it is likely that ArgOfPericenter and MeanAnomaly are both incorrect\n"
				<< "\tas a result. You will need to fix them. If an object's Argument of Pericenter is undefined in US, then MeanAnomaly here is likely incorrect.\n"
				<< "\tThere are other cases where the numbers look fine in US but are still broken here. And I'm sure there are even more cases that I couldn't test for.\n\n"

				<< "\tWhat I'm saying is, you need to check the AscendingNode, ArgOfPericenter, and MeanAnomaly for any object that looks wrong. To be honest, even\n"
				<< "\tif it looks fine in Space Engine, it could still be wrong in this file. The more complicated your system is, the more inaccurate the numbers will\n"
				<< "\tbe here. I hate that it has to be this way, but there is nothing I can do about it.\t*/\n\n\n";
			PrintFile(planetFile, *root);
			planetFile.close();

			object.clear();
		}
		closedir(dir);
	}
	else
	{
		// could not open directory
		std::cout << "\nThere was an error opening the \"input\" folder! Make sure it exists!\n";
		return 0;
	}

	std::cout << "\n Conversion complete! Look in the \"output\" folder\n to find your files.\n ";
	return 0;
}

void BuildHierarchy(std::list<Object>& binary, std::list<Object>::iterator& b)
{
	int size = object.size();
	double force, dist, max_f, min_d;
	Object *attractor = NULL, *neighbor = NULL;

	// sort by descending mass
	std::sort(object.begin(), object.end(), [](Object const& one, Object const& two){ return ( one.mass > two.mass ); } );

	if (size>1)
	{
		for (int i = 0; i < size; i++)
		{
			max_f = 0.0;
			min_d = 1.79769e308;
			for (int j = 0; j < size; j++)
			{
				if (i==j) continue;
				// identify strongest attractor
				force = object.at(j).mass / pow(Distance(object.at(i), object.at(j)), 2); // G, m constants
				if (max_f < force && (i==0 || object.at(i).mass < object.at(j).mass))
				{
					max_f = force;
					attractor = &object.at(j);
				}
				// identify closest neighbor
				dist = Distance(object.at(i), object.at(j));
				if (dist < min_d)
				{
					min_d = dist;
					neighbor = &object.at(j);
				}
			}
			if (attractor == neighbor){ // else { guaranteed: attractor > neighbor }

				object.at(i).mass > attractor->mass ? Bond(binary, b, object.at(i), *attractor) : Bond(binary, b, *attractor, object.at(i));
			}
			else if (object.at(i).mass > neighbor->mass && object.at(i).mass > attractor->mass)
			{
				Bond(binary, b, object.at(i), *attractor);
				Distance(*neighbor, object.at(i)) < Distance(*neighbor, *attractor) ? Bond(binary, b, object.at(i), *neighbor) : Bond(binary, b, *attractor, *neighbor);
			}
			else if (object.at(i).mass < neighbor->mass && object.at(i).mass < attractor->mass)
			{
				Bond(binary, b, *attractor, *neighbor);

				neighbor->hillSphereRadius = Distance(*neighbor, *attractor) * (cbrt(neighbor->mass / (3 * attractor->mass)));
				dist = Distance(object.at(i), *neighbor);
				if (dist < neighbor->hillSphereRadius)
					Bond(binary, b, *neighbor, object.at(i));
				else // force leftovers to pick a dominant object
					Bond(binary, b, *attractor, object.at(i));
			}
			else // neighbor < obj < attractor
			{
				Bond(binary, b, *attractor, object.at(i));

				object.at(i).hillSphereRadius = Distance(object.at(i), *attractor) * (cbrt(object.at(i).mass / (3 * attractor->mass)));
				dist = Distance(object.at(i), *neighbor);
				if (dist < object.at(i).hillSphereRadius)
					Bond(binary, b, object.at(i), *neighbor);
			}
		}
	}

	// as long as all objects are connected up the hierarchy, this will find the root
	root = &object.at(0);
	while (root->parent != NULL)
		root = root->parent;
	root->partner = root->parent = root;
}

void Bond(std::list<Object>& binary, std::list<Object>::iterator& b, Object& dom, Object& sub)
{
	if (sub.parent == &dom || &dom == sub.partner)
		return;

	// whether bonding happens with the massive body or with its barycenter (if partner != parent it's bound to ITS OWN bary)
	if ((dom.parent != NULL && dom.parent->type == "Barycenter") && dom.parent != dom.partner)
	{
		double forceObj = (sub.mass * dom.mass) / pow(Distance(sub, dom), 2);
		double forceBary = (sub.mass * dom.parent->mass) / pow(Distance(sub, *dom.parent), 2);
		if (forceBary > forceObj)
		{
			if (Distance(sub, *dom.parent) > Distance(*dom.partner, *dom.parent))
			{
				Bond(binary, b, *dom.parent, sub);
				return;
			}
		}
	}

	// same as above for when Sub already has a Bary with something else than Dom
	if ((sub.parent != NULL && sub.parent->type == "Barycenter") && sub.parent != sub.partner)
	{
		double forceObj = (dom.mass * sub.mass) / pow(Distance(dom, sub), 2);
		double forceBary = (dom.mass * sub.parent->mass) / pow(Distance(dom, *sub.parent), 2);
		if (forceBary > forceObj)
		{
			Bond(binary, b, dom, *sub.parent);
			return;
		}
	}

	// whether bonding is direct or requires a barycenter
	if ((dom.type == "Star" && sub.type == "Star") || (sub.mass / dom.mass) > 0.01)
	{
		// create barycenter with heavier object first
		Object temp;
		temp.type = "Barycenter";
		temp.parent = dom.parent;
		temp.partner = dom.partner;
		temp.child.push_back(&dom);
		temp.child.push_back(&sub);
		binary.push_back(temp);
		b++;

		if(dom.partner != NULL && dom.partner != dom.parent)
			dom.partner->partner = &*b;
		if (dom.parent != NULL || sub.parent != NULL)
		{
			for (int i=0; i < b->parent->child.size(); i++)
			{
				if (b->parent->child.at(i) == &dom || b->parent->child.at(i) == &sub)
				b->parent->child.erase(b->parent->child.begin() + i--);
			}
			b->parent->child.push_back(&*b);
		}

		dom.partner = &sub;
		sub.partner = &dom;
		dom.parent = sub.parent = &*b;
		UpdateBinary(*b, dom, sub);
	}
	else
	{
		sub.partner = sub.parent = &dom;
		sub.parent->child.push_back(&sub);
	}
}

// binaries of binaries have to be corrected when hierarchy changes
void UpdateBinary(Object& bin, Object& A, Object& B)
{
	bin.mass = A.mass + B.mass;
	bin.name = A.mass > B.mass ? (A.name + "-" + B.name) : (B.name + "-" + A.name);
	// a barycenter is "stellar" if one of its components is a star
	bin.isStar = A.isStar || B.isStar;
	// This essentially finds the average of the weighted vectors to determine the position of the barycenter
	double AposRatio = A.mass / bin.mass,
		BposRatio = B.mass / bin.mass;
	StateVect AposScaled = Scale(AposRatio, A.position),
		BposScaled = Scale(BposRatio, B.position);
	bin.position = Add(AposScaled, BposScaled);

	// this also averages the velocity vectors, weighted based on distance from barycenter, to find the velocity of the barycenter
	double distTotal = Distance(A, B),
		AvelRatio = (Distance(B, bin) / distTotal),
		BvelRatio = (Distance(A, bin) / distTotal);
	StateVect AvelScaled = Scale(AvelRatio, A.velocity),
		BvelScaled = Scale(BvelRatio, B.velocity);
	bin.velocity = Add(AvelScaled, BvelScaled);

	// should bin.tilt and bin.angularVelocity be calculated for use in orbits of children?

	if (bin.parent != NULL && bin.parent->type == "Barycenter" && bin.partner->partner == &bin)
		UpdateBinary(*bin.parent, bin, *bin.partner);
}

void CalcOrbit(Object& obj)
{
	/* The function starts with the bottom-most object of the hierarchy and works its way to the top. 
		mu is set manually for each object depending on the parent. Then, the object position and velocity
		subtract that of the parent's to center the object around it's parent. Then it uses those values with
		the functions from the paper. The only objects this does not work for are binary objects. */
	for (int i = 0; i < obj.child.size(); i++)
		CalcOrbit(*obj.child.at(i));

	if (&obj == root)
		return;

	double mu = G * obj.parent->mass;

	obj.position = Subtract(obj.position, obj.parent->position);
	obj.velocity = Subtract(obj.velocity, obj.parent->velocity);

	StateVect momentVect;
	// calculate the momentum vector h
	momentVect = CrossProduct(obj.position, obj.velocity);

	// calculate inclination
	obj.inclination = (acos(momentVect.z / momentVect.magnitude()));

	// calculate eccentricity vector and eccentricity
	StateVect eccentVect, VcrossH = CrossProduct(obj.velocity, momentVect);
	eccentVect.x = ((VcrossH.x / mu) - (obj.position.x / obj.position.magnitude()));
	eccentVect.y = ((VcrossH.y / mu) - (obj.position.y / obj.position.magnitude()));
	eccentVect.z = ((VcrossH.z / mu) - (obj.position.z / obj.position.magnitude()));
	obj.eccentricity = eccentVect.magnitude();

	// calculate true anomaly
	double Tanomaly;
	if (DotProduct(obj.position, obj.velocity) >= 0.0)
		Tanomaly = (acos(DotProduct(eccentVect, obj.position) / (eccentVect.magnitude() * obj.position.magnitude())));
	else
		Tanomaly = ((2 * PI) - acos(DotProduct(eccentVect, obj.position) / (eccentVect.magnitude() * obj.position.magnitude())));

	// calculate vector n / for ascending node
	StateVect K(0, 0, 1), n = CrossProduct(K, momentVect);

	// calculate long of Ascending Node
	if (n.y >= 0.0)
		obj.longOfAscNode = (acos(n.x / n.magnitude()) + PI); // Why I have to add Pi here and subtract it from the next equation is unknown to me... But it seems to work...?
	else
		obj.longOfAscNode = ( (2 * PI) - acos(n.x / n.magnitude()) - PI);

	// calculate argOfPeriapsis
	if (eccentVect.z >= 0.0)
		obj.argOfPeriapsis = (acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));
	else
		obj.argOfPeriapsis = ((2 * PI) - acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));

	// calculate eccentric anomaly
	double E = (2 * atan((tan(Tanomaly / 2)) / (sqrt((1 + obj.eccentricity) / (1 - obj.eccentricity)))));

	// calculate mean anomaly
	obj.meanAnomaly = (E - (obj.eccentricity * sin(E)));

	/*
		Obliquity doesn't work and I don't know how to fix it

		If YOU know how to fix it, send me a message and we will figure it out. 
		Here is the problem:

		Orientation from the simulation file is represented as four numbers, I assume
		to be a quaternion of the form r, i, j, k. It looks something like this:
		"Orientation":"-0.03764766;-0.8111074;-0.03670029;0.5811343",

		Those four numbers need to be turned into obliquity relative to the
		orbital plane somehow. I was using this pseudo-code from the US team
		to get an idea:

			up = orbitnormal
			rotationaxis = normalize(AngularVelocity)
			degree = 180 / PI
			world = normalize(rotationaxis * Orientation)
			obliquity = (PI - acos(dot(up, world))) * degree

		The problem is, I don't know if I am reading it correctly as I have
		no clue what to do for Orientation being a quaternion. Does
		rotationaxis scale orientation? Does orientation need to be
		in euler angles first? idk and honestly idc anymore

		If you can read that pseudo-code and know how to slove this, 
		please tell me and we will get obliquity working!
	*/
	// calculates obliquity
	StateVect orbitNormal = Normalize(momentVect);
	StateVect rotationAxis = Normalize(obj.angularVelocity);
	StateVect orientation = To_EulerAngles(obj.orientation);
	StateVect world = Scale(rotationAxis.magnitude(), orientation);
	world = Normalize(world);
	obj.obliquity = (PI - acos(DotProduct(orbitNormal, world)));
	obj.obliquity = To_Degree(obj.obliquity);

	//StateVect tiltVect = RotateVector(obj.orientation, obj.angularVelocity);
	//obj.obliquity = acos(DotProduct(tiltVect, momentVect) / (tiltVect.magnitude() * momentVect.magnitude()));
	//obj.obliquity = abs(180 - To_Degree(obj.obliquity));

	return;
}

void CalcMoreOrbit(Object& obj)
{
	for (int i = 0; i < obj.child.size(); i++)
		CalcMoreOrbit(*obj.child.at(i));

	if (&obj == root)
		return;

	double mu, mu2 = 0, semimajor2 = 0;
	mu = G * obj.partner->mass;

	// calculate semi-major axis
	obj.semimajor = ( 1 / ((2 / obj.position.magnitude()) - (pow(obj.velocity.magnitude(), 2) / mu)) );
	// calculate partner's semi-major axis if they're binaries
	if (&obj == obj.partner->partner)
	{
		mu2 = G * obj.mass;
		semimajor2 = ( 1 / ((2 / obj.partner->position.magnitude()) - (pow(obj.partner->velocity.magnitude(), 2) / mu2)) );
	}

	obj.period = ( (2 * PI) * (sqrt( pow(obj.semimajor + semimajor2, 3) / ( mu + mu2 ) )) );
	obj.period /= 3.1556952e7; // converts sec to years;

	obj.semimajor /= 1.496e11; // converts m to AU
	//obj.semimajor /= 1000; // converts m to km

	// hydrostatic equilibrium is determined by shape, but SU2 doesn't seem to store this information
	// large asteroids that pass the mass test are still not round because of frequent collisions (their orbits aren't clear)
	obj.isRound = (obj.mass > HYDROUPPER) || (obj.mass > HYDROLOWER && ClearanceCheck(obj));

	return;
}

bool ClearanceCheck(Object& obj)
{
	// Pi discriminant for orbital clearance (a mass range of 10e23 can be used instead,
	// which would make 4 solar system moons dwarf planets if they weren't moons)
	if (obj.partner != obj.parent) // orbits own barycenter
		return ( K_KGKGAU * (obj.mass / (pow(obj.parent->parent->mass, 5.0/2.0) * pow(obj.parent->semimajor, 9.0/8.0))) > 1 );
	else
		return ( K_KGKGAU * (obj.mass / (pow(obj.parent->mass, 5.0/2.0) * pow(obj.semimajor, 9.0/8.0))) > 1 );
}

StateVect To_EulerAngles(Quaternion q) {
	StateVect angles;

	// roll (x-axis rotation)
	double sinr_cosp = 2 * (q.w * q.x + q.y * q.z);
	double cosr_cosp = 1 - 2 * (q.x * q.x + q.y * q.y);
	angles.x = atan2(sinr_cosp, cosr_cosp);

	// pitch (y-axis rotation)
	double sinp = 2 * (q.w * q.y - q.z * q.x);
	if (abs(sinp) >= 1)
		angles.y = copysign(PI / 2, sinp); // use 90 degrees if out of range
	else
		angles.y = asin(sinp);

	// yaw (z-axis rotation)
	double siny_cosp = 2 * (q.w * q.z + q.x * q.y);
	double cosy_cosp = 1 - 2 * (q.y * q.y + q.z * q.z);
	angles.z = atan2(siny_cosp, cosy_cosp);

	return angles;
}

// rotate vector based on a quaternion's rotation matrix
StateVect RotateVector(Quaternion& q, StateVect& vect)
{
	StateVect result;

	result.x = vect.x * (q.w * q.w + q.x * q.x - q.y * q.y - q.z * q.z);
	result.x += vect.y * (2 * q.x * q.y - 2 * q.w * q.z);
	result.x += vect.z * (2 * q.x * q.z + 2 * q.w * q.y);

	result.y = vect.x * (2 * q.x * q.y + 2 * q.w * q.z);
	result.y += vect.y * (q.w * q.w - q.x * q.x + q.y * q.y - q.z * q.z);
	result.y += vect.z * (2 * q.y * q.z - 2 * q.w * q.x);

	result.z = vect.x * (2 * q.x * q.z - 2 * q.w * q.y);
	result.z += vect.y * (2 * q.y * q.z + 2 * q.w * q.x);
	result.z += vect.z * (q.w * q.w - q.x * q.x - q.y * q.y + q.z * q.z);

	return result;
}

double Distance(Object& A, Object& B)
{
	// finds the distance between two objects in 3d space
	double x, y, z;
	x = pow(B.position.x - A.position.x, 2);
	y = pow(B.position.y - A.position.y, 2);
	z = pow(B.position.z - A.position.z, 2);
	return sqrt(x + y + z);
}

StateVect CrossProduct(StateVect& A, StateVect& B)
{
	StateVect cross;
	cross.x = ((A.y * B.z) - (B.y * A.z));
	cross.y = ((A.z * B.x) - (B.z * A.x));
	cross.z = ((A.x * B.y) - (B.x * A.y));
	return cross;
}

double DotProduct(StateVect& A, StateVect& B)
{
	return ((A.x * B.x) + (A.y * B.y) + (A.z * B.z));
}

StateVect Scale(double scaler, StateVect& A)
{
	StateVect temp;
	temp.x = A.x * scaler;
	temp.y = A.y * scaler;
	temp.z = A.z * scaler;
	return temp;
}

StateVect Subtract(StateVect& A, StateVect& B)
{
	StateVect temp;
	temp.x = A.x - B.x;
	temp.y = A.y - B.y;
	temp.z = A.z - B.z;
	return temp;
}

StateVect Add(StateVect& A, StateVect& B)
{
	StateVect temp;
	temp.x = A.x + B.x;
	temp.y = A.y + B.y;
	temp.z = A.z + B.z;
	return temp;
}

StateVect Normalize(StateVect& a)
{
	double mag = a.magnitude();
	StateVect temp;
	temp.x = a.x / mag;
	temp.y = a.y / mag;
	temp.z = a.z / mag;
	return temp;
}

double To_Degree(double& number)
{
	return (number * (180 / PI));
}

void GetData(std::ifstream& inputFile)
{
	std::string holder;

	while (true)
	{
		Object temp;
		temp.parent = NULL;
		temp.partner = NULL;
		temp.isStar = false;
		temp.ironMass = 0;
		temp.waterMass = 0;
		temp.hydrogenMass = 0;
		double compositionMassTotal = 0.0;

		// finds name
		while (inputFile >> holder && holder != "\"$type\":\"Body\",");
		if (holder == "}")
			return;
		std::getline(inputFile, holder);
		std::getline(inputFile, holder);
		holder.erase(0, 8);
		holder.erase(holder.size() - 2, 2);
		temp.name = holder;

		// find albedo
		while (inputFile >> holder && !(holder.find("\"Albedo\":") + 1) && !(holder.find("\"Id\":") + 1));
		if (holder.find("\"Id\":") + 1) // If an object does not have an albedo, it skips the next few items because they don't exist either
		{
			temp.temp = 0;
			temp.albedo = 0;
			temp.greenhouse = 0;
			temp.surfacePressure = 0;
			temp.luminosity = 0;
			temp.roughness = 0;
			temp.isStar = false;
			goto NoAlbedo;
		}
		holder.erase(0, 9);
		temp.albedo = std::stod(holder, &sz);

		// find temperature
		while (inputFile >> holder && !(holder.find("\"SurfaceTemperature\":") + 1));
		holder.erase(0, 21);
		temp.temp = std::stoi(holder, &sz);

		// find greenhouse
		while (inputFile >> holder && !(holder.find("\"GreenhouseEffect\":") + 1));
		holder.erase(0, 19);
		temp.greenhouse = std::stod(holder, &sz);

		// find surface pressure
		while (inputFile >> holder && !(holder.find("\"SurfacePressure\":") + 1));
		holder.erase(0, 18);
		temp.surfacePressure = std::stod(holder, &sz);
		temp.surfacePressure /= 101325; // converts pa to atm

		// find luminosity
		while (inputFile >> holder && !(holder.find("\"Luminosity\":") + 1));
		holder.erase(0, 13);
		temp.luminosity = std::stod(holder, &sz);
		temp.luminosity /= 3.827e26; // converts watts to solar lum

		// find color stuff
		while (inputFile >> holder && !(holder.find("\"Colors\":[") + 1));
		if (inputFile >> holder && holder != "],") // object uses low/mid/high colors
		{
			holder.erase(0, 6);
			temp.colorUpland.r = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorUpland.g = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorUpland.b = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorUpland.a = std::stod(holder, &sz);

			holder.erase(0, 14);
			temp.colorLowland.r = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorLowland.g = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorLowland.b = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorLowland.a = std::stod(holder, &sz);

			holder.erase(0, 14);
			temp.colorBeach.r = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorBeach.g = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorBeach.b = std::stod(holder, &sz);
			inputFile >> holder;
			temp.colorBeach.a = std::stod(holder, &sz);	
		}
		else // object has no colors or uses color layers
		{
			while (inputFile >> holder && !(holder.find("\"Colors\":[") + 1));
			if (inputFile >> holder && holder != "],") // object uses color layers
			{
				holder.erase(0, 6);
				temp.colorLayer0.r = std::stod(holder, &sz);
				inputFile >> holder;
				temp.colorLayer0.g = std::stod(holder, &sz);
				inputFile >> holder;
				temp.colorLayer0.b = std::stod(holder, &sz);
				inputFile >> holder;
				temp.colorLayer0.a = std::stod(holder, &sz);
				holder.erase(0, 7);

				if (holder != "],")
				{
					holder.erase(0, 7);
					temp.colorLayer1.r = std::stod(holder, &sz);
					inputFile >> holder;
					temp.colorLayer1.g = std::stod(holder, &sz);
					inputFile >> holder;
					temp.colorLayer1.b = std::stod(holder, &sz);
					inputFile >> holder;
					temp.colorLayer1.a = std::stod(holder, &sz);
					holder.erase(0, 7);

					if (holder != "],")
					{
						holder.erase(0, 7);
						temp.colorLayer2.r = std::stod(holder, &sz);
						inputFile >> holder;
						temp.colorLayer2.g = std::stod(holder, &sz);
						inputFile >> holder;
						temp.colorLayer2.b = std::stod(holder, &sz);
						inputFile >> holder;
						temp.colorLayer2.a = std::stod(holder, &sz);
						holder.erase(0, 7);

						if (holder != "],")
						{
							holder.erase(0, 7);
							temp.colorLayer3.r = std::stod(holder, &sz);
							inputFile >> holder;
							temp.colorLayer3.g = std::stod(holder, &sz);
							inputFile >> holder;
							temp.colorLayer3.b = std::stod(holder, &sz);
							inputFile >> holder;
							temp.colorLayer3.a = std::stod(holder, &sz);
							holder.erase(0, 7);

							if (holder != "],")
							{
								holder.erase(0, 7);
								temp.colorLayer4.r = std::stod(holder, &sz);
								inputFile >> holder;
								temp.colorLayer4.g = std::stod(holder, &sz);
								inputFile >> holder;
								temp.colorLayer4.b = std::stod(holder, &sz);
								inputFile >> holder;
								temp.colorLayer4.a = std::stod(holder, &sz);
								holder.erase(0, 7);

								if (holder != "],")
								{
									holder.erase(0, 7);
									temp.colorLayer5.r = std::stod(holder, &sz);
									inputFile >> holder;
									temp.colorLayer5.g = std::stod(holder, &sz);
									inputFile >> holder;
									temp.colorLayer5.b = std::stod(holder, &sz);
									inputFile >> holder;
									temp.colorLayer5.a = std::stod(holder, &sz);
									holder.erase(0, 7);

									if (holder != "],")
									{
										holder.erase(0, 7);
										temp.colorLayer6.r = std::stod(holder, &sz);
										inputFile >> holder;
										temp.colorLayer6.g = std::stod(holder, &sz);
										inputFile >> holder;
										temp.colorLayer6.b = std::stod(holder, &sz);
										inputFile >> holder;
										temp.colorLayer6.a = std::stod(holder, &sz);
										holder.erase(0, 7);

										if (holder != "],")
										{
											holder.erase(0, 7);
											temp.colorLayer7.r = std::stod(holder, &sz);
											inputFile >> holder;
											temp.colorLayer7.g = std::stod(holder, &sz);
											inputFile >> holder;
											temp.colorLayer7.b = std::stod(holder, &sz);
											inputFile >> holder;
											temp.colorLayer7.a = std::stod(holder, &sz);
										}
									}
								}
							}
						}
					}
				}
			}
		}

		// Find mass composition, if it exists
		while (inputFile >> holder &&
			holder != "\"Iron\":{" && holder != "\"Water\":{" && holder != "\"Hydrogen\":{" &&
			(!(holder.find("\"ElevationToRadiusRatio\":") + 1)));
		// if the mass is 100% silicate, the next part is skipped over
		if (holder == "\"Iron\":{")
		{
			// find iron mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.ironMass = std::stod(holder, &sz);

			while (inputFile >> holder &&
				holder != "\"Water\":{" && holder != "\"Hydrogen\":{" &&
				(!(holder.find("\"ElevationToRadiusRatio\":") + 1)));
		}
		if (holder == "\"Water\":{")
		{
			// find water mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.waterMass = std::stod(holder, &sz);

			while (inputFile >> holder &&
				holder != "\"Hydrogen\":{" &&
				(!(holder.find("\"ElevationToRadiusRatio\":") + 1)));
		}
		if (holder == "\"Hydrogen\":{")
		{
			// find hydrogen mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.hydrogenMass = std::stod(holder, &sz);

			while (inputFile >> holder && (!(holder.find("\"ElevationToRadiusRatio\":") + 1)));
		}

		// find surface roughness
		holder.erase(0, 25);
		temp.roughness = std::stod(holder, &sz);

	NoAlbedo:;

		// find mass, isStar, type, shape
		while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
		holder.erase(0, 7);
		temp.mass = std::stod(holder, &sz);
		//temp.mass = (temp.mass / (5.9736 * pow(10, 24))); // converts kg to earth masses
		if (temp.mass > STARCUTOFF)
		{
			temp.isStar = true;
			temp.type = "Star";
		}
		else
		{
			temp.isStar = false;
			temp.type = ""; // will be dealt with after hierarchy is found
		}

		// converts composition into %'s to be printed for SE
		compositionMassTotal = temp.hydrogenMass + temp.waterMass + temp.ironMass;
		temp.silicateMass = temp.mass - compositionMassTotal;
		temp.hydrogenMass /= temp.mass / 100;
		temp.silicateMass /= temp.mass / 100;
		temp.waterMass /= temp.mass / 100;
		temp.ironMass /= temp.mass / 100;

		// find radius
		while (inputFile >> holder && !(holder.find("\"Radius\":") + 1));
		holder.erase(0, 9);
		temp.radius = std::stod(holder, &sz);
		temp.radius /= 1000; // m to km
		temp.roughness *= temp.radius; // get the elevation span

		std::string x, y, z;

		// find orientation and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Orientation\":") + 1));
		holder.erase(0, 15);
		x = holder.substr(holder.find(";"));
		x.erase(0, 1);
		y = x.substr(x.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.orientation.w = std::stod(holder, &sz); // quaternion element r
		temp.orientation.x = std::stod(x, &sz); // quaternion element i
		temp.orientation.y = std::stod(y, &sz); // quaternion element j
		temp.orientation.z = std::stod(z, &sz); // quaternion element k

		// finds angular velocity
		while (inputFile >> holder && !(holder.find("\"AngularVelocity\":") + 1));
		holder.erase(0, 19);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.angularVelocity.x = std::stod(holder, &sz); // x
		temp.angularVelocity.y = std::stod(y, &sz); // y
		temp.angularVelocity.z = std::stod(z, &sz); // z
		temp.rotationPeriod = temp.angularVelocity.magnitude();
		temp.rotationPeriod = ((2 * PI) / temp.rotationPeriod); // calculates rot period
		temp.rotationPeriod /= 3600; // sec to hours

		// find position and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Position\":") + 1));
		holder.erase(0, 12);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.position.x = std::stod(holder, &sz); // x
		temp.position.z = std::stod(y, &sz); // y
		temp.position.y = std::stod(z, &sz); // z

		// find velocity and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Velocity\":") + 1));
		holder.erase(0, 12);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.velocity.x = std::stod(holder, &sz); // x
		temp.velocity.z = std::stod(y, &sz); // y
		temp.velocity.y = std::stod(z, &sz); // z

		// find category
		while (inputFile >> holder && !(holder.find("\"Category\":") + 1));
		holder.erase(0, 12);
		holder.erase(holder.size() - 1, 1);
		temp.class_ = holder;

		// add object to global vector, ignoring US2 barycenters (which have R=1m, M=1kg) and nameless fragments
		if (temp.name != "" && temp.radius != 0.001)
			object.push_back(temp);
	}
}

// evaluate type according to orbital relationships and mass (as defined by IAU / used by SE)
void Typifier(Object& obj)
{
	if (&obj == root || (obj.partner != obj.parent && (obj.parent == root && obj.mass > obj.partner->mass)))
	{
		obj.type = "Star"; // even an asteroid-size planemo is considered a "star" by SE
		return;
	}

	if (obj.isStar)
		return;
	else if (obj.partner->isStar)
	{
		if (ClearanceCheck(obj))
			obj.type = "Planet";
		else if (obj.isRound)
			obj.type = "DwarfPlanet";
		else
			// mean comet eccentricity is 0.7, while most e=0.4 SSSO's are comets or asteroids of comet origins
			obj.type = (obj.eccentricity > 0.4) ? "Comet" : "Asteroid";
	}
	else if (obj.partner->isRound) // there's no separate type for Moon-moons in SE (yet)
	{
		if (obj.isRound)
			if (obj.mass < obj.partner->mass)
				obj.type = "Moon";
			else
				obj.type = ClearanceCheck(obj) ? "Planet" : "DwarfPlanet";
		else
			obj.type = "DwarfMoon";
	}
	else // incl. binary comet or asteroid
		obj.type = (obj.eccentricity > 0.4) ? "Comet" : "Asteroid";
}

void Classifier(Object& obj)
{
	if (obj.type == "Star")
	{
		double density = 0.75 * obj.mass / (PI * pow(1000 * obj.radius, 3));
		if (obj.class_ == "blackhole")
			obj.class_ = "X";
		else if (density > 1e16)
			obj.class_ = "Q";
		else if (density > 1e8)
			obj.class_ = "WD";
		else
			obj.class_ = "";
	}
	else if (obj.isRound) // boundaries currently used by SE, see their blog
	{
		if (obj.hydrogenMass > 0.01)
		{
			if (obj.hydrogenMass < 25)
				obj.class_ = "Neptune";
			else
				obj.class_ = "Jupiter";
		}
		else
		{
			if (obj.ironMass > 50)
				obj.class_ = "Ferria";
			else if (obj.waterMass > 25)
				obj.class_ = "Aquaria";
			else
				obj.class_ = "Terra";
		}
	}
	else
		obj.class_ = "Asteroid";
}

void PrintFile(std::ofstream& f, Object & o)
{
	f << o.type << "\t\t\t\t\t\"" << o.name << "\""
		<< "\n{"
		<< "\n\tParentBody\t\t\t\"" << o.parent->name << "\"";
	if (&o == root && o.type == "Barycenter")
	{
		f << "\n}\n\n";
		for (int i = 0; i < o.child.size(); i++)
			PrintFile(f, *o.child.at(i));
		return;
	}
	else if (&o == root && o.type == "Star")
	{
		f << "\n\tClass\t\t\t\t\"" << o.class_ << "\""
			<< "\n\tLum\t\t\t\t\t" << o.luminosity
			<< "\n\tTeff\t\t\t\t" << o.temp
			<< "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
			<< "\n\tRadius\t\t\t\t" << o.radius
			<< "\n\tRotationPeriod\t\t" << o.rotationPeriod
			<< "\n\tObliquity\t\t\t" << ((o.obliquity < -180) ? 0.0 : o.obliquity)
			<< "\n}\n\n";
		for (int i = 0; i < o.child.size(); i++)
			PrintFile(f, *o.child.at(i));
		return;
	}

	if (o.type != "Barycenter")
	{
		if (o.class_ != "")
			f << "\n\tClass\t\t\t\t\"" << o.class_ << "\"";
		f << "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
			<< "\n\tRadius\t\t\t\t" << o.radius
			<< "\n\tRotationPeriod\t\t" << o.rotationPeriod
			<< "\n\tObliquity\t\t\t" << o.obliquity;
		if (o.type == "Star")
			f << "\n\tLum\t\t\t\t\t" << o.luminosity
				<< "\n\tTeff\t\t\t\t" << o.temp;
		else
		{
			f << "\n\tAlbedoBond\t\t\t" << o.albedo
				<< "\n\n\tComposition"
				<< "\n\t{"
				<< "\n\t\tHydrogen\t" << o.hydrogenMass
				<< "\n\t\tHelium\t\t0"
				<< "\n\t\tSilicates\t" << o.silicateMass
				<< "\n\t\tCarbides\t0"					// helium/carbide output added for the sake of the user
				<< "\n\t\tIces\t\t" << o.waterMass
				<< "\n\t\tMetals\t\t" << o.ironMass
				<< "\n\t}";

			f << "\n\n\tSurface"
				<< "\n\t{";
			if (o.class_ == "Jupiter" || o.class_ == "Neptune")
			{
				f << "\n\t\tcolorLayer0\t\t" << o.colorLayer0.Print();
				if (o.colorLayer1.r >= 0)
				{
					f << "\n\t\tcolorLayer1\t\t" << o.colorLayer1.Print();
					if (o.colorLayer2.r >= 0)
					{
						f << "\n\t\tcolorLayer2\t\t" << o.colorLayer2.Print();
						if (o.colorLayer3.r >= 0)
						{
							f << "\n\t\tcolorLayer3\t\t" << o.colorLayer3.Print();
							if (o.colorLayer4.r >= 0)
							{
								f << "\n\t\tcolorLayer4\t\t" << o.colorLayer4.Print();
								if (o.colorLayer5.r >= 0)
								{
									f << "\n\t\tcolorLayer5\t\t" << o.colorLayer5.Print();
									if (o.colorLayer6.r >= 0)
									{
										f << "\n\t\tcolorLayer6\t\t" << o.colorLayer6.Print();
										if (o.colorLayer7.r >= 0)
										{
											f << "\n\t\tcolorLayer7\t\t" << o.colorLayer7.Print();
										}
									}
								}
							}
						}
					}
				}
			}
			else
				f << "\n\t\tcolorShelf\t\t" << o.colorBeach.Print()
				<< "\n\t\tcolorBeach\t\t" << o.colorBeach.Print()
				<< "\n\t\tcolorDesert\t\t" << o.colorLowland.Print()
				<< "\n\t\tcolorLowland\t" << o.colorLowland.Print()
				<< "\n\t\tcolorUpland\t\t" << o.colorUpland.Print()
				<< "\n\t\tcolorRock\t\t" << o.colorUpland.Print();
			f << "\n\t}";
		}
			
	}

	if (o.type != "Barycenter" && o.hydrogenMass < 0.01)
	{
		f << "\n\n\tAtmosphere"
			<< "\n\t{"
			<< "\n\t\tPressure\t\t" << o.surfacePressure
			<< "\n\t\tGreenhouse\t\t" << o.greenhouse
			<< "\n\t}";
	}

	f << "\n\n\tOrbit"
		<< "\n\t{"
		<< "\n\t\tRefPlane\t\t\"Equator\""
		<< "\n\t\tSemiMajorAxis\t" << o.semimajor
		<< "\n\t\tPeriod\t\t\t" << o.period
		<< "\n\t\tEccentricity\t" << o.eccentricity
		<< "\n\t\tInclination\t\t" << To_Degree(o.inclination)
		<< "\n\t\tAscendingNode\t" << To_Degree(o.longOfAscNode)
		<< "\n\t\tArgOfPericenter\t" << To_Degree(o.argOfPeriapsis)
		<< "\n\t\tMeanAnomaly\t\t" << To_Degree(o.meanAnomaly)
		<< "\n\t}";

	f << "\n}\n\n";

	for (int i = 0; i < o.child.size(); i++)
		PrintFile(f, *o.child.at(i));
	return;
}
