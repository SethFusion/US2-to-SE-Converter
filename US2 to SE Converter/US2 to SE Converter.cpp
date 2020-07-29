#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>

#include <algorithm> // for sort()
#include <limits> // for max double

/*
	Listen here...

	Before you start judging this code for how ugly it is, just remember that I'm no physics major and
	this shit really took it out of me. Yeah, there are probably some optimizations this garbage could
	use. In fact I know there are parts of this I only wrote for testing, and then left in because it
	actually worked. I'll even say this: it is the ugliest program I've ever written. But you know what?
	It gets the job done, for the general case, and that's all I need it to do.
*/

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
		return (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
	}
};

struct Object
{
	std::string name, type, class_;
	Object* parent;
    Object* partner; // for binaries and objects of similar mass, otherwise it's == parent
	std::vector<Object*> child;
	double mass, radius;

	// orbit stuff
	StateVect position, velocity, angularVelocity, orientation;
	double semimajor, period, inclination, eccentricity, argOfPeriapsis, longOfAscNode, meanAnomaly, hillSphereRadius;
	double rotationPeriod, obliquity;

	// general
	int temp;
	double ironMass, waterMass, hydrogenMass, surfacePressure, greenhouse;

	// used for stars
	bool isStar;
	double luminosity;

	// asteroid?
	bool isAsteroid;
};
std::string::size_type sz;
std::vector<Object> object;
Object* root;

double G; // gravitation constant
const double PI = 3.1415926535; // it's pi you idiot

void GetData(std::ifstream&);
void PrintFile(std::ofstream&, Object&);
void Bond(std::list<Object>&, std::list<Object>::iterator&, Object*, Object*);
void CreateBinary(std::list<Object>&, Object&, Object&, Object*);

// math functions
double Distance(Object&, Object&);
void CalcOrbit(Object&);
void CalcMoreOrbit(Object& obj);
StateVect CrossProduct(StateVect&, StateVect&);
double DotProduct(StateVect&, StateVect&);
StateVect Scale(double, StateVect&);
double Determinate(std::vector<double>&);
StateVect Subtract(StateVect&, StateVect&);
StateVect Add(StateVect&, StateVect&);
StateVect Normalize(StateVect&);

int main()
{
	std::string systemName, starFileName, planetFileName, holder;

	std::ifstream inputFile("input/simulation.json");
	if (!inputFile)
	{
		std::cout << "\nThere was an error opening the simulation file! Make sure it is placed in the \"input\" folder!\n";
		return 0;
	}

	// finds the "Name": part of the simulation file and uses it for star and planet .sc files
	while (inputFile >> holder && holder != "},");
	while (inputFile >> holder && holder != "},");
	std::getline(inputFile, holder);
	std::getline(inputFile, holder);
	starFileName = holder;
	starFileName.erase(0, 8);
	starFileName.erase(starFileName.size() - 2, 2);
	planetFileName = systemName = starFileName;
	starFileName = "output/" + starFileName + " Star.sc";
	planetFileName = "output/" + planetFileName + " Planet.sc";

	// finds gravitational constant
	while (inputFile >> holder && !(holder.find("\"Gravity\":") + 1));
	holder.erase(0, 10);
	G = std::stod(holder, &sz);

	// inputs every object into the global vector
	GetData(inputFile);
	inputFile.close();

    Object *attractor = NULL, *neighbor = NULL;
    double force, dist, max_f, min_d;

	std::list<Object> binary;
	Object systemO;
	systemO.name = systemName + " System";
	binary.push_back(systemO);
	std::list<Object>::iterator b = binary.begin(); // binary object counter

	int size = object.size();

	std::sort(object.begin(), object.end(), [](Object const& one, Object const& two){ return ( one.mass > two.mass ); } ); // supposed to be pretty quick in c++?

    for (int i = 0; i < size; i++)
    {
std::cout << "  " << object.at(i).name << "\n";
    }

    for (int i = 0; i < size-2; i++)
    {
        max_f = 0.0;
        min_d = std::numeric_limits<double>::max();
        for (int j = 0; j < size-2; j++)
        {
            if (i==j) continue;
            // identify strongest attractor
            force = (object.at(i).mass * object.at(j).mass) / pow(Distance(object.at(i), object.at(j)), 2);
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
std::cout << "1" << "\n";
            object.at(i).mass > attractor->mass ? Bond(binary, b, &object.at(i), attractor) : Bond(binary, b, attractor, &object.at(i));
        }
        else if (object.at(i).mass > neighbor->mass && object.at(i).mass > attractor->mass)
        {
std::cout << "2" << "\n";
            Bond(binary, b, &object.at(i), attractor);
            Distance(*neighbor, object.at(i)) < Distance(*neighbor, *attractor) ? Bond(binary, b, &object.at(i), neighbor) : Bond(binary, b, attractor, neighbor);
        }
        else if (object.at(i).mass < neighbor->mass && object.at(i).mass < attractor->mass)
        {
std::cout << "3" << "\n";
            Bond(binary, b, attractor, neighbor);

            neighbor->hillSphereRadius = Distance(*neighbor, *attractor) * (cbrt(neighbor->mass / (3 * attractor->mass)));
            dist = Distance(object.at(i), *neighbor);
            if (dist < neighbor->hillSphereRadius)
                Bond(binary, b, neighbor, &object.at(i));
            else // force leftovers to pick a dominant object
                Bond(binary, b, attractor, &object.at(i));
        }
        else // neighbor < obj < attractor
        {
std::cout << "4" << "\n";
            Bond(binary, b, attractor, &object.at(i));

            object.at(i).hillSphereRadius = Distance(object.at(i), *attractor) * (cbrt(object.at(i).mass / (3 * attractor->mass)));
            dist = Distance(object.at(i), *neighbor);
            if (dist < object.at(i).hillSphereRadius)
                Bond(binary, b, &object.at(i), neighbor);
        }
    }

	// as long as all objects are connected up the hierarchy, this will find the root
	root = &object.at(0);
	while (root->parent != NULL)
		root = root->parent;

    for (int i = 0; i < size; i++)
    {
        if (object.at(i).parent == NULL)
        {
            object.at(i).parent = root;
std::cout << "5.   " << object.at(i).name << "\n";
        }
        if (object.at(i).partner == NULL)
{
            object.at(i).partner = root;
std::cout << "6.   " << object.at(i).name << "\n";}
    }
std::cout << "7: end.\n";

    // this is where the magic happens
	CalcOrbit(*root);
	// requires partners (binaries) to already have some orbit info
	CalcMoreOrbit(*root);

	std::ofstream starFile(starFileName.c_str());
	starFile << "StarBarycenter\t\t\"" << binary.begin()->name << "\"\n{}\n";
	starFile.close();

	root->parent = &*binary.begin();
	std::ofstream planetFile(planetFileName.c_str());
	PrintFile(planetFile, *root);

	std::cout << "\n Conversion complete! Look in the \"output\" folder\n to find your files.\n ";

	planetFile.close();
	return 0;
}

void Bond(std::list<Object>& binary, std::list<Object>::iterator& b, Object* dom, Object* sub)
{
    if (sub->parent == dom || dom == sub->partner)
        return;

    // whether bonding happens with the massive body or with its barycenter (if partner != parent it's bound to ITS OWN bary)
    if ((dom->parent != NULL && dom->parent->type == "Barycenter") && dom->parent != dom->partner)
    {
        double forceObj = (sub->mass * dom->mass) / pow(Distance(*sub, *dom), 2);
        double forceBary = (sub->mass * dom->parent->mass) / pow(Distance(*sub, *dom->parent), 2);
        if (forceBary > forceObj)
            dom = dom->parent;
    }

    // whether bonding is direct or requires a barycenter
    if ((dom->type == "Star" && sub->type == "Star") || (sub->mass / dom->mass) > 0.01)
    {
        CreateBinary(binary, *dom, *sub, dom->partner);
        dom->partner = sub;
        sub->partner = dom;
        b++;

        if (dom->parent != NULL)
        {
            b->parent = dom->parent;
            for (int i=0; i < b->parent->child.size(); i++)
            {
                if (b->parent->child.at(i) == dom || b->parent->child.at(i) == sub)
                b->parent->child.erase(b->parent->child.begin() + i--);
            }
            b->parent->child.push_back(&*b);
        }
        dom->parent = sub->parent = &*b;
std::cout << "b.   " << dom->name << "+" << sub->name << "    partners:   " << dom->partner->name << sub->partner->name << "\n";
    }
    else
    {
        sub->partner = sub->parent = dom;
        sub->parent->child.push_back(sub);
std::cout << "c.   " << dom->name << "+" << sub->name << "\n";
    }
}

void CreateBinary(std::list<Object>& binaryList, Object& A, Object& B, Object* C)
{
	Object temp;
	temp.isAsteroid = false;
	temp.isStar = false;
	temp.name = (A.name + "-" + B.name);
	temp.type = "Barycenter";
std::cout << "B.   " << temp.name << "\n";
	temp.mass = A.mass + B.mass;
	temp.parent = A.parent;
    if (C != NULL && A.parent != NULL)
        temp.partner = temp.parent->type == "Barycenter" ? C : temp.parent;
	// This essentially finds the average of the weighted vectors to determine the position of the barycenter
	double AposRatio = A.mass / temp.mass,
		BposRatio = B.mass / temp.mass;
	StateVect AposScaled = Scale(AposRatio, A.position),
		BposScaled = Scale(BposRatio, B.position);
	temp.position = Add(AposScaled, BposScaled);

	// this also averages the velocity vectors, weighted based on distance from barycenter, to find the velocity of the barycenter
	double distTotal = Distance(A, B),
		AvelRatio = (Distance(B, temp) / distTotal),
		BvelRatio = (Distance(A, temp) / distTotal);
	StateVect AvelScaled = Scale(AvelRatio, A.velocity),
		BvelScaled = Scale(BvelRatio, B.velocity);
	temp.velocity = Add(AvelScaled, BvelScaled);

	temp.child.push_back(&A);
	temp.child.push_back(&B);
	//temp.mass = 5.75e18;
	binaryList.push_back(temp);
}

void PrintFile(std::ofstream& f, Object & o)
{

	f << o.type << "\t\t\t\t\t\"" << o.name << "\""
		<< "\n{";
	if (&o == root && o.type == "Barycenter")
	{
		f << "\n\tParentBody\t\t\t\"" << o.parent->name << "\""
			<< "\n}\n\n";
		for (int i = 0; i < o.child.size(); i++)
			PrintFile(f, *o.child.at(i));
		return;
	}
	else if (&o == root && o.type == "Star")
	{
		f << "\n\tParentBody\t\t\t\"" << o.parent->name << "\""
			<< "\n\tLum\t\t\t\t\t" << o.luminosity
			<< "\n\tTeff\t\t\t\t" << o.temp
			<< "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
			<< "\n\tRadius\t\t\t\t" << o.radius
			<< "\n\tRotationPeriod:\t\t" << o.rotationPeriod
            //<< "\n\tObliquity:\t\t" << o.obliquity
			<< "\n}\n\n";
		for (int i = 0; i < o.child.size(); i++)
			PrintFile(f, *o.child.at(i));
		return;
	}
	else
		f << "\n\tParentBody\t\t\t\"" << o.parent->name << "\"";

	if (o.type != "Barycenter")
		f << "\n\tClass\t\t\t\t\"" << o.class_ << "\""
		<< "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
		<< "\n\tRadius\t\t\t\t" << o.radius
		<< "\n\tRotationPeriod:\t\t" << o.rotationPeriod;

	//f << "\n\tObliquity:\t\t" << o.obliquity
	f << "\n\n\tOrbit"
		<< "\n\t{"
		<< "\n\t\tRefPlane\t\t\"Equator\""
		<< "\n\t\tSemiMajorAxis\t" << o.semimajor
		<< "\n\t\tPeriod\t\t\t" << o.period
		<< "\n\t\tEccentricity\t" << o.eccentricity
		<< "\n\t\tInclination\t\t" << o.inclination
		<< "\n\t\tAscendingNode\t" << o.longOfAscNode
		<< "\n\t\tArgOfPericenter\t" << o.argOfPeriapsis
		<< "\n\t\tMeanAnomaly\t\t" << o.meanAnomaly
		<< "\n\t}";

	if (o.type != "Barycenter")
		f << "\n\n\tAtmosphere"
		<< "\n\t{"
		<< "\n\t\tPressure\t\t" << o.surfacePressure
		<< "\n\t\tGreenhouse\t\t" << o.greenhouse
		<< "\n\t}";

	f << "\n}\n\n";

	for (int i = 0; i < o.child.size(); i++)
		PrintFile(f, *o.child.at(i));
	return;
}

void CalcOrbit(Object& obj)
{
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
	obj.inclination = (obj.inclination * (180 / PI));
	obj.inclination -= 90; // convert to equator

	// calculate eccentricity vector and eccentricity
	StateVect eccentVect, VcrossH = CrossProduct(obj.velocity, momentVect);
	eccentVect.x = ((VcrossH.x / mu) - (obj.position.x / obj.position.magnitude()));
	eccentVect.y = ((VcrossH.y / mu) - (obj.position.y / obj.position.magnitude()));
	eccentVect.z = ((VcrossH.z / mu) - (obj.position.z / obj.position.magnitude()));
	obj.eccentricity = eccentVect.magnitude();

	// calculate vector n / for ascending node
	StateVect one(0.0, 0.0, 1.0), two(momentVect.y * -1.0, momentVect.x, 0.0);
	StateVect n = CrossProduct(one, two);

	// calculate true anomaly
	double Tanomaly;
	if (DotProduct(obj.position, obj.velocity) >= 0.0)
		Tanomaly = (acos(DotProduct(eccentVect, obj.position) / (eccentVect.magnitude() * obj.position.magnitude())));
	else
		Tanomaly = ((2 * PI) - acos(DotProduct(eccentVect, obj.position) / (eccentVect.magnitude() * obj.position.magnitude())));

	// calculate eccentric anomaly
	double E = (2 * atan( (tan(Tanomaly / 2)) / (sqrt((1 + obj.eccentricity) / (1 - obj.eccentricity)) )));

	// calculate long of Ascending Node
	if (n.y >= 0.0)
		obj.longOfAscNode = (acos(n.x / n.magnitude()));
	else
		obj.longOfAscNode = ( (2 * PI) - acos(n.x / n.magnitude()));
	obj.longOfAscNode = (obj.longOfAscNode * (180 / PI));
	obj.longOfAscNode -= 90; // convert to equator

	// calculate argOfPeriapsis
	if (eccentVect.z >= 0.0)
		obj.argOfPeriapsis = (acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));
	else
		obj.argOfPeriapsis = ((2 * PI) - acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));

	// Calculates Obliquity
	obj.obliquity = (acos(DotProduct(obj.angularVelocity, obj.parent->angularVelocity) / (obj.parent->angularVelocity.magnitude()*obj.parent->angularVelocity.magnitude()) ));
	obj.obliquity = (obj.obliquity * (180 / PI)) - obj.inclination;
	//obj.obliquity -= 90; // convert to equator


/*
	StateVect left, forward, parentTemp;
	parentTemp.x = (obj.argOfPeriapsis - obj.parent->position.x);
	parentTemp.y = (obj.argOfPeriapsis - obj.parent->position.y);
	parentTemp.z = (obj.argOfPeriapsis - obj.parent->position.z);
	left = Normalize(parentTemp);
	forward = CrossProduct(left, momentVect);

	StateVect rotationAxis, world;
	rotationAxis = Normalize(obj.angularVelocity);
	// world =
	obj.obliquity = (PI - acos(DotProduct(momentVect, world)));
	// argument of obliquity

*/



	obj.argOfPeriapsis = (obj.argOfPeriapsis * (180 / PI)); // convert to degree
	obj.argOfPeriapsis -= 90; // convert to equator

	// calculate mean anomaly
	obj.meanAnomaly = (E - (obj.eccentricity * sin(E)));
	obj.meanAnomaly = (obj.meanAnomaly * (180 / PI));
	obj.meanAnomaly -= 90; // convert to equator

	return;
}

void CalcMoreOrbit(Object& obj)
{
	for (int i = 0; i < obj.child.size(); i++)
		CalcMoreOrbit(*obj.child.at(i));

	if (&obj == root)
		return;

	double mu, mu2, semimajor2;
	mu = G * obj.partner->mass;
	mu2 = G * obj.mass;

	// calculate semi-major axis
	obj.semimajor = ( 1 / ((2 / obj.position.magnitude()) - (pow(obj.velocity.magnitude(), 2) / mu)) );
	// calculate partner's semi-major axis (even if very small)
	semimajor2 = ( 1 / ((2 / obj.partner->position.magnitude()) - (pow(obj.partner->velocity.magnitude(), 2) / mu2)) );

	obj.period = ( (2 * PI) * (sqrt( pow(obj.semimajor + semimajor2, 3) / ( mu + mu2 ) )) );
	obj.period /= 3.1556952e7; // converts sec to years;

	obj.semimajor /= 1.496e11; // converts m to AU
	//obj.semimajor /= 1000; // converts m to km

	return;
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
	std::vector<double> i, j, k;

	i.push_back(A.y);
	i.push_back(A.z);
	i.push_back(B.y);
	i.push_back(B.z);

	j.push_back(A.x);
	j.push_back(A.z);
	j.push_back(B.x);
	j.push_back(B.z);

	k.push_back(A.x);
	k.push_back(A.y);
	k.push_back(B.x);
	k.push_back(B.y);

	cross.x = Determinate(i);
	cross.y = (Determinate(j) * -1.0);
	cross.z = Determinate(k);
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

double Determinate(std::vector<double>& vect)
{
	return ((vect.at(0) * vect.at(3)) - (vect.at(1) * vect.at(2)));
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






void GetData(std::ifstream& inputFile)
{
	std::string holder;

	while (true)
	{
		Object temp;
		temp.parent = NULL;
		temp.isAsteroid = false;
		temp.isStar = false;
		temp.ironMass = 0;
		temp.waterMass = 0;
		temp.hydrogenMass = 0;

		// finds name
		while (inputFile >> holder && holder != "\"$type\":\"Body\",");
		if (holder == "}")
			return;
		std::getline(inputFile, holder);
		std::getline(inputFile, holder);
		holder.erase(0, 8);
		holder.erase(holder.size() - 2, 2);
		temp.name = holder;

		// find temperature
		while (inputFile >> holder && !(holder.find("\"SurfaceTemperature\":") + 1) && !(holder.find("\"Id\":") + 1));
		if (holder.find("\"Id\":") + 1) // If an object does not have a temperature, it skips the next few items because they don't exist either
		{
			temp.isAsteroid = true;
			temp.temp = 0;
			temp.greenhouse = 0;
			temp.surfacePressure = 0;
			temp.luminosity = 0;
			temp.isStar = false;
			goto NoTemperature;
		}
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
		temp.luminosity = (temp.luminosity / (3.827 * pow(10, 26))); // converts watts to solar lum

		// determine if it is a star
		while (inputFile >> holder && !(holder.find("\"StarType\":") + 1));
		holder.erase(0, 11);
		temp.isStar = std::stoi(holder, &sz);

		// find molten level
		while (inputFile >> holder && !(holder.find("\"MoltenLevel\":") + 1));
		// Find mass composition, if it exists
		while (inputFile >> holder &&
			holder != "\"Iron\":{" && holder != "\"Water\":{" && holder != "\"Hydrogen\":{" &&
			(!(holder.find("\"Id\":") + 1)));
		// if the mass is 100% silicate, the next part is skipped over
		if (holder == "\"Iron\":{")
		{
			// find iron mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.ironMass = std::stod(holder, &sz);
			//temp.ironMass = (temp.ironMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

			while (inputFile >> holder &&
				holder != "\"Water\":{" && holder != "\"Hydrogen\":{" &&
				(!(holder.find("\"Id\":") + 1)));
		}
		if (holder == "\"Water\":{")
		{
			// find water mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.waterMass = std::stod(holder, &sz);
			//temp.waterMass = (temp.waterMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

			while (inputFile >> holder &&
				holder != "\"Hydrogen\":{" &&
				(!(holder.find("\"Id\":") + 1)));
		}
		if (holder == "\"Hydrogen\":{")
		{
			// find hydrogen mass
			while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
			holder.erase(0, 7);
			temp.hydrogenMass = std::stod(holder, &sz);
			//temp.hydrogenMass = (temp.hydrogenMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

			while (inputFile >> holder && (!(holder.find("\"Id\":") + 1)));
		}
	NoTemperature:;

		// find mass
		while (inputFile >> holder && !(holder.find("\"Mass\":") + 1));
		holder.erase(0, 7);
		temp.mass = std::stod(holder, &sz);
		//temp.mass = (temp.mass / (5.9736 * pow(10, 24))); // converts kg to earth masses

		// find radius
		while (inputFile >> holder && !(holder.find("\"Radius\":") + 1));
		holder.erase(0, 9);
		temp.radius = std::stod(holder, &sz);
		temp.radius /= 1000; // m to km

		std::string y, z;
		/*
		// find orientation and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Orientation\":") + 1));
		holder.erase(0, 15);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.orientation.x = std::stod(holder, &sz); // x
		temp.orientation.y = std::stod(y, &sz); // y
		temp.orientation.z = std::stod(z, &sz); // z
		// still one value left
		*/

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
		temp.rotationPeriod = ((2 * PI) / abs(temp.rotationPeriod)); // calculates rot period
		temp.rotationPeriod /= 3600; // sec to hours

		// find position and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Position\":") + 1));
		holder.erase(0, 12);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.position.x = std::stod(holder, &sz); // x
		temp.position.y = std::stod(y, &sz); // y
		temp.position.z = std::stod(z, &sz); // z

		// find velocity and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Velocity\":") + 1));
		holder.erase(0, 12);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.velocity.x = std::stod(holder, &sz); // x
		temp.velocity.y = std::stod(y, &sz); // y
		temp.velocity.z = std::stod(z, &sz); // z

		// find category / type
		while (inputFile >> holder && !(holder.find("\"Category\":") + 1));
		holder.erase(0, 12);
		holder.erase(holder.size() - 1, 1);
		temp.type = holder;

		// uses category to determine what kind of object this is
		if (temp.isStar == true || temp.type == "star")
		{
			temp.type = "Star";
			temp.class_ = "";
		}
		else if (temp.type == "planet")
		{
			temp.type = "Planet";
			if (temp.hydrogenMass / temp.mass > 0.01)
			{
				if (temp.waterMass / temp.mass > 0.3)
					temp.class_ = "Neptune";
				else
					temp.class_ = "Jupiter";
			}
			else
			{
				if (temp.waterMass / temp.mass > 0.05)
					temp.class_ = "Aquaria";
				else if (temp.ironMass / temp.mass > 0.3)
					temp.class_ = "Ferria";
				else
					temp.class_ = "Terra";
			}
		}
		else if (temp.type == "moon")
		{
			if (temp.radius > 200)
				temp.type = "Moon";
			else
				temp.type = "DwarfMoon";

			if (temp.waterMass / temp.mass > 0.05)
				temp.class_ = "Aquaria";
			else if (temp.ironMass / temp.mass > 0.3)
				temp.class_ = "Ferria";
			else
				temp.class_ = "Terra";
		}
		else if (temp.type == "sso")
		{
			if (temp.isAsteroid || temp.radius < 300)
			{
				temp.type = "Asteroid";
				temp.class_ = "Asteroid";
			}
			else
			{
				temp.type = "DwarfPlanet";

				if (temp.waterMass / temp.mass > 0.05)
					temp.class_ = "Aquaria";
				else if (temp.ironMass / temp.mass > 0.3)
					temp.class_ = "Ferria";
				else
					temp.class_ = "Terra";
			}
		}
		else if (temp.type == "blackhole")
		{
			temp.type = "Star";
			temp.class_ = "X";
		}

		// adds object to global vector
		object.push_back(temp);
	}
}

