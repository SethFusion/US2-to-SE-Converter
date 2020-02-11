
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

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
	double mass, radius;
	int temp;
	double age, ironMass, waterMass, hydrogenMass, surfacePressure, greenhouse;

	// orbit stuff
	Object* parent;
	std::vector<Object*> child;
	StateVect position, velocity;
	double semimajor, inclination, eccentricity, argOfPeriapsis, longOfAscNode, meanAnomaly, hillSphereRadius;

	// used for stars
	bool isStar;
	double luminosity;

	// asteroid?
	bool isAsteroid;
};
std::string::size_type sz;
std::vector<Object> object;
double G; // gravitation constant
const double PI = 3.1415926535; // it's pi you idiot

void GetData(std::ifstream&);
void PrintObject(std::ofstream&, Object&);
void CreateBinary(std::vector<Object>&, Object&, Object&);

// math functions
double Distance(Object&, Object&);
void CalcOrbit(Object&);
StateVect CrossProduct(StateVect&, StateVect&);
double DotProduct(StateVect&, StateVect&);
double Determinate(std::vector<double>&);
StateVect Subtract(StateVect&, StateVect&);

int main()
{
	std::string systemName, starFileName, planetFileName, holder;
	bool choice = true; // true = full system, false = closed system

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

	GetData(inputFile);
	inputFile.close();


	std::vector<Object> star, planet, moon, binary;
	int size = object.size();
	for (int i = 0; i < size; i++)
	{
		if (object.at(i).type == "Moon" || object.at(i).type == "DwarfMoon" || object.at(i).type == "Asteroid")
			moon.push_back(object.at(i));
		else if (object.at(i).type == "Planet" || object.at(i).type == "DwarfPlanet")
			planet.push_back(object.at(i));
		else if (object.at(i).type == "Star")
			star.push_back(object.at(i));
	}

	for (int i = 0; i < planet.size(); i++)
	{
		double F = 0.0; // force of the current star on the current planet
		int parent = 0; // positon of the most attractive star
		for (int j = 0; j < star.size(); j++)
		{
			double f = (G * pow(planet.at(i).mass, 2) * pow(star.at(j).mass, 2)) / pow(Distance(planet.at(i), star.at(j)), 2);
			// checks the current star's attractive force with the largest known force
			// If the current star's force is larger than any others, the parent
			// is set to the current star
			if (f > F)
			{
				F = f;
				parent = j;
			}
		}
		planet.at(i).parent = &star.at(parent);
		star.at(parent).child.push_back(&planet.at(i));
		planet.at(i).hillSphereRadius = Distance(planet.at(i), star.at(parent)) * (cbrt(planet.at(i).mass / (3 * star.at(parent).mass)));
	}

	// looks for binary planet systems
	for (int i = 0; i < planet.size(); i++)
	{
		for (int j = 0; j < planet.size(); j++)
		{
			double dist = Distance(planet.at(i), planet.at(j));
			if (dist < planet.at(i).hillSphereRadius && dist != 0)
			{
				planet.at(i).parent = &planet.at(j);
				planet.at(j).parent = &planet.at(i);
			}
		}
	}

	// so bassically this uses the psudo-hill sphere radius of the planet to check if any moons are within that radius
	// If they are within it, it likely means that planet is the moons parentbody
	for (int i = 0; i < moon.size(); i++)
	{
		int parent = 0;
		for (int j = 0; j < planet.size(); j++)
		{
			if (Distance(moon.at(i), planet.at(j)) < planet.at(j).hillSphereRadius)
			{
				parent = j;
				moon.at(i).parent = &planet.at(j);
			}
		}
		moon.at(i).hillSphereRadius = Distance(moon.at(i), planet.at(parent)) * (cbrt(moon.at(i).mass / (3 * planet.at(parent).mass)));
		// If the moon fails to find a parent body, it will search the star list
		// to for whatever star is pulling on it the most
		if (moon.at(i).parent == NULL)
		{
			double F = 0.0;
			for (int j = 0; j < star.size(); j++)
			{
				double f = (G * pow(moon.at(i).mass, 2) * pow(star.at(j).mass, 2)) / pow(Distance(moon.at(i), star.at(j)), 2);
				if (f > F)
				{
					F = f;
					parent = j;
				}
			}
			moon.at(i).parent = &star.at(parent);
		}
	}

	// this is where the magic happens
	for (int i = 0; i < moon.size(); i++)
	{
		CalcOrbit(moon.at(i));
		if (moon.at(i).type == "Asteroid" && moon.at(i).eccentricity > 0.5)
			moon.at(i).type = "Comet";
	}
	for (int i = 0; i < planet.size(); i++)
	{
		CalcOrbit(planet.at(i));
	}
	


	std::cout << " Do you want an empty star file for your system?"
		<< "\n\n"
		<< " Type 1 for Yes or 0 for No, then press enter.\n\n ";
	std::cin >> choice;

	if (choice)
	{
		std::ofstream starFile(starFileName.c_str());

		starFile << "StarBarycenter\t\t\"" << systemName << " System\"\n{}\n";
		starFile.close();
	}

	std::ofstream planetFile(planetFileName.c_str());
	for (int i = 0; i < star.size(); i++)
	{
		PrintObject(planetFile, star.at(i));
	}
	for (int i = 0; i < planet.size(); i++)
	{
		PrintObject(planetFile, planet.at(i));
	}
	for (int i = 0; i < moon.size(); i++)
	{
		PrintObject(planetFile, moon.at(i));
	}

	std::cout << "\n Conversion complete! Look in the \"output\" folder\n to find your files.\n ";
	std::cin >> choice;

	planetFile.close();
	return 0;
}

void PrintObject(std::ofstream& f, Object & o)
{
	
	if (!(o.type == "Star"))
	{
		f << o.type << "\t\t\t\t\t\"" << o.name << "\""
			<< "\n{"
			<< "\n\tParentBody\t\t\t\"" << o.parent->name << "\""
			<< "\n\tClass\t\t\t\t\"" << o.class_ << "\""
			<< "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
			<< "\n\tRadius\t\t\t\t" << o.radius
			//<< "\n\tAge\t\t\t\t\t" << o.age
			<< "\n\n\tOrbit"
			<< "\n\t{"
			<< "\n\t\tRefPlane\t\t\"Ecliptic\""
			<< "\n\t\tSemiMajorAxis\t" << o.semimajor
			<< "\n\t\tEccentricity\t" << o.eccentricity
			<< "\n\t\tInclination\t\t" << o.inclination
			<< "\n\t\tAscendingNode\t" << o.longOfAscNode
			<< "\n\t\tArgOfPericenter\t" << o.argOfPeriapsis
			<< "\n\t\tMeanAnomaly\t\t" << o.meanAnomaly
			<< "\n\t}"
			<< "\n\n\tAtmosphere"
			<< "\n\t{"
			<< "\n\t\tPressure\t\t" << o.surfacePressure
			<< "\n\t\tGreenhouse\t\t" << o.greenhouse
			<< "\n\t}"
			<< "\n}\n\n";
		return;
	}

	f << o.type << "\t\t\t\t\t\"" << o.name << "\""
		<< "\n{"
		<< "\n\tParentBody\t\t\t\"aEarth\""
		<< "\n\tLum\t\t\t\t\t" << o.luminosity
		<< "\n\tMass\t\t\t\t" << o.mass / (5.9736 * pow(10, 24))
		<< "\n\tRadius\t\t\t\t" << o.radius
		<< "\n\tTeff\t\t\t\t" << o.temp
		<< "\n\tAge\t\t\t\t\t" << o.age
		<< "\n}\n\n";
	return;
}

void CreateBinary(std::vector<Object>& binaryVect, Object& A, Object& B)
{
	Object temp;
	temp.name = (A.name + "-" + B.name);
	temp.type = "Barycenter";
	temp.mass = A.mass + B.mass;
	temp.parent = A.parent;

	double a = Distance(A, B);
	

	temp.child.push_back(&A);
	temp.child.push_back(&B);
}

void CalcOrbit(Object& obj)
{
	double mu;
	mu = G * obj.parent->mass;

	obj.position = Subtract(obj.position, obj.parent->position);
	obj.velocity = Subtract(obj.velocity, obj.parent->velocity);


	StateVect momentVect;
	// caculate the momentum vector h
	momentVect = CrossProduct(obj.position, obj.velocity);

	// calculate inclination
	obj.inclination = (acos(momentVect.z / momentVect.magnitude()));
	obj.inclination = (obj.inclination * (180 / PI));

	// calcuate eccentricity vector and eccentricity
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

	// calculate argOfPeriapsis
	if (eccentVect.z >= 0.0)
		obj.argOfPeriapsis = (acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));
	else
		obj.argOfPeriapsis = ((2 * PI) - acos(DotProduct(n, eccentVect) / (n.magnitude() * eccentVect.magnitude())));
	obj.argOfPeriapsis = (obj.argOfPeriapsis * (180 / PI));

	// calculate mean anomaly
	obj.meanAnomaly = (E - (obj.eccentricity * sin(E)));
	obj.meanAnomaly = (obj.meanAnomaly * (180 / PI));

	// calulate semi-major axis
	obj.semimajor = ( 1 / ((2 / obj.position.magnitude()) - (pow(obj.velocity.magnitude(), 2) / mu)) );
	obj.semimajor /= 1.496e11;
	//obj.semimajor /= 1000;

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

		// find age
		while (inputFile >> holder && !(holder.find("\"Age\":") + 1));
		holder.erase(0, 6);
		temp.age = std::stod(holder, &sz);
		temp.age /= 31557600000000000; // converst seconds to gigayears

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
