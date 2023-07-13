#include <testutilscxx.h>

template<typename T>
void add_all(std::vector<T> &vector1, std::vector<T> vector2)
{
	vector1.insert(vector1.end(), vector2.begin(), vector2.end());
}

void
assert_equals(CuTest *tc,
              const TubularMesh &t1,
              const TubularMesh &t2)
{
	CuAssertIntEquals(tc, t1.tubularSegments(), t2.tubularSegments());
	CuAssertIntEquals(tc, t1.radialSegments(), t2.radialSegments());
	CuAssertDblEquals(tc, t1.radius(), t2.radius(), POINT_EPSILON);

	CuAssertIntEquals(tc, 0, t1.vertices().size() % 3);
	CuAssertIntEquals(tc, 0, t1.normals().size() % 3);
	CuAssertIntEquals(tc, 0, t1.tangents().size() % 4);
	CuAssertIntEquals(tc, 0, t1.uvs().size() % 2);

	CuAssertIntEquals(tc, t1.vertices().size(), t2.vertices().size());
	CuAssertIntEquals(tc, t1.normals().size(), t2.normals().size());
	CuAssertIntEquals(tc, t1.tangents().size(), t2.tangents().size());
	CuAssertIntEquals(tc, t1.uvs().size(), t2.uvs().size());

	size_t vSize, iSize;
	std::vector<real> t1_vectors = t1.vertices();
	add_all(t1_vectors, t1.normals());
	add_all(t1_vectors, t1.tangents());
	add_all(t1_vectors, t1.uvs());
	std::vector<real> t2_vectors = t2.vertices();
	add_all(t2_vectors, t2.normals());
	add_all(t2_vectors, t2.tangents());
	add_all(t2_vectors, t2.uvs());
	TubularMesh::size(t1.tubularSegments(), t1.radialSegments(), vSize, iSize);
	for (int i = 0; i < vSize; i++)
		CuAssertDblEquals(tc, t1_vectors.at(i), t2_vectors.at(i), POINT_EPSILON);

	std::vector<size_t> t1_indices = t1.indices();
	std::vector<size_t> t2_indices = t2.indices();
	for (int i = 0; i < iSize; i++)
		CuAssertIntEquals(tc, t1_indices.at(i), t2_indices.at(i));
}

void
tubularmesh_copy_ctor(CuTest *tc)
{
	// Given
	BSpline spline(5, 3);
	spline.setControlPoints({
			100.0, 200.0, 0.0,   /* P1 */
			150.0, 220.0, 10.0,  /* P2 */
			190.0, 120.0, 50.0,  /* P3 */
			260.0, 70.0,  30.0,  /* P4 */
			300.0, 200.0, 20.0 /* P5 */
	});
	TubularMesh mesh = spline.computeRMF(spline.uniformKnotSeq(51)).createTube(50, 8, 10.0);

	// When
	TubularMesh copy(mesh);

	// Then
	assert_equals(tc, mesh, copy);
}

void
tubularmesh_move_ctor(CuTest *tc)
{
	// Given
	BSpline spline(5, 3);
	spline.setControlPoints({
			100.0, 200.0, 0.0,   /* P1 */
			150.0, 220.0, 10.0,  /* P2 */
			190.0, 120.0, 50.0,  /* P3 */
			260.0, 70.0,  30.0,  /* P4 */
			300.0, 200.0, 20.0 /* P5 */
	});
	TubularMesh mesh;
	mesh = spline.computeRMF(spline.uniformKnotSeq(51)).createTube(50, 8, 10.0);
	TubularMesh toBeMoved(mesh);

	// When
	TubularMesh move(std::move(toBeMoved));

	// Then
	assert_equals(tc, mesh, move);
	CuAssertIntEquals(tc, 0, toBeMoved.tubularSegments());
	CuAssertIntEquals(tc, 0, toBeMoved.radialSegments());
	CuAssertDblEquals(tc, 0, toBeMoved.radius(), POINT_EPSILON);
}

void
tubularmesh_copy_assign(CuTest *tc)
{
	// Given
	BSpline spline(5, 3);
	spline.setControlPoints({
			100.0, 200.0, 0.0,   /* P1 */
			150.0, 220.0, 10.0,  /* P2 */
			190.0, 120.0, 50.0,  /* P3 */
			260.0, 70.0,  30.0,  /* P4 */
			300.0, 200.0, 20.0 /* P5 */
	});
	TubularMesh mesh = spline.computeRMF(spline.uniformKnotSeq(51)).createTube(50, 8, 10.0);
	TubularMesh copy;

	// When
	copy = mesh;

	// Then
	assert_equals(tc, mesh, copy);
}

void
tubularmesh_move_assign(CuTest *tc)
{
	// Given
	BSpline spline(5, 3);
	spline.setControlPoints({
		100.0, 200.0, 0.0,   /* P1 */
		150.0, 220.0, 10.0,  /* P2 */
		190.0, 120.0, 50.0,  /* P3 */
		260.0, 70.0,  30.0,  /* P4 */
		300.0, 200.0, 20.0 /* P5 */
	});
	TubularMesh mesh = spline.computeRMF(spline.uniformKnotSeq(51)).createTube(50, 8, 10.0);
	TubularMesh toBeMoved(mesh);
	TubularMesh move;

	// When
	move = std::move(toBeMoved);

	// Then
	assert_equals(tc, mesh, move);
	CuAssertIntEquals(tc, 0, toBeMoved.tubularSegments());
	CuAssertIntEquals(tc, 0, toBeMoved.radialSegments());
	CuAssertDblEquals(tc, 0, toBeMoved.radius(), POINT_EPSILON);
}

CuSuite *
get_tubularmesh_suite()
{
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, tubularmesh_copy_ctor);
	SUITE_ADD_TEST(suite, tubularmesh_move_ctor);
	SUITE_ADD_TEST(suite, tubularmesh_copy_assign);
	SUITE_ADD_TEST(suite, tubularmesh_move_assign);
	return suite;
}
