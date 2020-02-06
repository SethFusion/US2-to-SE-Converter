
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

struct Object
{
	std::string name, type, class_;
	double mass, radius;
	int temp;
	std::vector<double> position;
	std::vector<double> velocity;
	double age, ironMass, waterMass, hydrogenMass, surfacePressure, greenhouse;

	//used for stars
	bool isStar;
	double luminosity;
};
std::vector<Object> object;

void getData(std::ifstream &);
void printObject(std::ofstream &, Object &);

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

	// This finds the "Name": part of the sinulation file and uses it for star and planet .sc files
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

	getData(inputFile);
	inputFile.close();

	std::cout << " Do you want an empty star file for your system?"
		<< "\n\n"
		<< " Type 1 for Yes or 0 for No, then press enter.\n\n ";
	std::cin >> choice;

	if (choice)
	{
		std::ofstream starFile(starFileName.c_str());

		starFile << "StarBarycenter\t\t\"" << systemName << "\"\n{}\n";
		starFile.close();
	}

	std::ofstream planetFile(planetFileName.c_str());
	for (int i = 0; i < object.size() - 1; i++)
	{
		printObject(planetFile, object.at(i));
	}

	std::cout << "\n Conversion complete! Look in the \"output\" folder\n to find your files.\n ";
	std::cin >> choice;

	planetFile.close();
	return 0;
}

void printObject(std::ofstream& f, Object & o)
{
	
	if (!(o.type == "Star"))
	{
		f << o.type << "\t\t\t\t\t\"" << o.name << "\""
			<< "\n{"
			<< "\n\tClass\t\t\t\t\"" << o.class_ << "\""
			<< "\n\tMass\t\t\t\t" << o.mass
			<< "\n\tRadius\t\t\t\t" << o.radius
			<< "\n\tAge\t\t\t\t\t" << o.age
			<< "\n\n\tAtmopshere"
			<< "\n\t{"
			<< "\n\t\tPressure\t\t" << o.surfacePressure
			<< "\n\t\tGreenhouse\t\t" << o.greenhouse
			<< "\n\t}"
			<< "\n}\n\n";
		return;
	}

	f << o.type << "\t\t\t\t\t\"" << o.name << "\""
		<< "\n{"
		<< "\n\tLum\t\t\t\t\t" << o.luminosity
		<< "\n\tMass\t\t\t\t" << o.mass
		<< "\n\tRadius\t\t\t\t" << o.radius
		<< "\n\tTeff\t\t\t\t" << o.temp
		<< "\n\tAge\t\t\t\t\t" << o.age
		<< "\n}\n\n";
	return;
}
void getData(std::ifstream& inputFile)
{
	std::string holder;
	std::string::size_type sz;

	while (true)
	{
		Object temp;
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
		if (holder.find("\"Id\":") + 1) // If an object does not have a temperature, it skips the next few iitems because they don't exist either
		{
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
			temp.ironMass = (temp.ironMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

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
			temp.waterMass = (temp.waterMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

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
			temp.hydrogenMass = (temp.hydrogenMass / (5.9736 * pow(10, 24))); // converts kg to earth masses

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
		temp.mass = (temp.mass / (5.9736 * pow(10, 24))); // converts kg to earth masses

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
		temp.position.push_back(std::stod(holder, &sz)); // x
		temp.position.push_back(std::stod(y, &sz)); // y
		temp.position.push_back(std::stod(z, &sz)); // z

		// find velocity and add x y z to vector
		while (inputFile >> holder && !(holder.find("\"Velocity\":") + 1));
		holder.erase(0, 12);
		y = holder.substr(holder.find(";"));
		y.erase(0, 1);
		z = y.substr(y.find(";"));
		z.erase(0, 1);
		temp.velocity.push_back(std::stod(holder, &sz)); // x
		temp.velocity.push_back(std::stod(y, &sz)); // y
		temp.velocity.push_back(std::stod(z, &sz)); // z

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
				if (temp.waterMass > 0.05)
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
			if (temp.radius > 200)
			{
				temp.type = "DwarfPlanet";

				if (temp.waterMass / temp.mass > 0.05)
					temp.class_ = "Aquaria";
				else if (temp.ironMass / temp.mass > 0.3)
					temp.class_ = "Ferria";
				else
					temp.class_ = "Terra";
			}
			else
			{
				temp.type = "Asteroid";
				temp.class_ = "Asteroid";
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
