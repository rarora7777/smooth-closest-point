#include <Eigen/Core>
#include <Phong.h>

int main()
{
	double vertices[] = {
		0, 1, 1, 0, 0, 1, 1, 0,
		0, 0, 1, 1, 1, 1, 0, 0,
		0, 0, 0, 0, 1, 1, 1, 1,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0
	};

	unsigned int triangles[] = {
		0, 0, 2, 2, 1, 1, 0, 0, 5, 5, 0, 0,
		2, 3, 3, 4, 2, 5, 7, 4, 4, 7, 6, 1,
		1, 2, 4, 5, 5, 6, 4, 3, 7, 6, 7, 6
	};

	auto phong = createPhongObject(vertices, 8, 8, triangles, 12);
	double point[] = { 1, 1.1, 1, 0, 0, 0, 0, 0 };
	float projection[4];
	bool res = project(phong, point, 0, projection);
	
	std::cout << int(projection[3]) << ' ' << projection[0] << ' ' << projection[1] << ' ' << projection[2] <<std::endl;

	deletePhongObject(phong);
	delete[] projection;
	return 0;
}