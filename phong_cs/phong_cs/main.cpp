#include <Eigen/Core>
#include <Phong.h>
#include <fstream>
#include <iostream>

int main()
{
	//double vertices[] = {
	//	0, 1, 1, 0, 0, 1, 1, 0,
	//	0, 0, 1, 1, 1, 1, 0, 0,
	//	0, 0, 0, 0, 1, 1, 1, 1,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0,
	//	0, 0, 0, 0, 0, 0, 0, 0
	//};

	//unsigned int triangles[] = {
	//	0, 0, 2, 2, 1, 1, 0, 0, 5, 5, 0, 0,
	//	2, 3, 3, 4, 2, 5, 7, 4, 4, 7, 6, 1,
	//	1, 2, 4, 5, 5, 6, 4, 3, 7, 6, 7, 6
	//};
	
	std::fstream file;
	file.open("../../../Assets/StreamingAssets/head_tri.txt");
	
	const int dim = 8;
	int nV, nF;
	file >> nV >> nF;
	
	std::vector<double> vertices(dim * nV);
	std::vector<unsigned int> triangles(3 * nF);
	
	for (int i = 0; i < nV; ++i)
	{
		double val;
		file >> val;
		file >> val;
		file >> val;
		
		for (int j = 0; j < dim; ++j)
		{
			file >> val;
			vertices[j * nV + i] = val;
		}
	}
            
	for (int i = 0; i < nF; ++i)
		for (int j = 0; j < 3; ++j)
		{
			unsigned int val;
			file >> val;
			triangles[j * nF + i] = val;
		}
	

	file.close();
	
	auto phong = createPhongObject(vertices.data(), nV, 8, triangles.data(), nF);

	std::cout << "Creation time: " << phong->initTime << std::endl;
	double point[] = { 1, 1.1, 1, 0, 0, 0, 0, 0 };
	float projection[4];
	bool res = project(phong, point, 0, projection);
	
	std::cout << int(projection[3]) << ' ' << projection[0] << ' ' << projection[1] << ' ' << projection[2] <<std::endl;

	deletePhongObject(phong);
	delete[] projection;
	return 0;
}