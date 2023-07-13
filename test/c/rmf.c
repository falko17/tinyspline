#include <testutils.h>

void
rmf_vectors_of_frames(CuTest *tc)
{
	___SETUP___
	tsBSpline spline = ts_bspline_init();
	tsBSpline derivative = ts_bspline_init();
	tsDeBoorNet net = ts_deboornet_init();
	tsReal knots[50], *result = NULL, dist, angle;
	tsFrame frames[50];
	size_t i;

	___GIVEN___
	C(ts_bspline_new_with_control_points(
		5, 3, 3, TS_CLAMPED, &spline, &status,
		100.0, 200.0, 0.0,   /* P1 */
		150.0, 220.0, 10.0,  /* P2 */
		190.0, 120.0, 50.0,  /* P3 */
		260.0, 70.0,  30.0,  /* P4 */
		300.0, 200.0, 20.0)) /* P5 */
	C(ts_bspline_derive(&spline,
	                    1,
	                    0,
	                    &derivative,
	                    &status))
	ts_bspline_uniform_knot_seq(&spline, 50, knots);

	___WHEN___
	C(ts_bspline_compute_rmf(&spline,
	                         knots,
	                         50,
	                         0,
	                         frames,
	                         &status))

	___THEN___
	for (i = 0; i < 50; i++) {
		/* Position */
		C(ts_bspline_eval(&spline, knots[i], &net, &status))
		C(ts_deboornet_result(&net, &result, &status))
		dist = ts_distance(result, frames[i].position, 3);
		CuAssertDblEquals(tc, 0, dist, POINT_EPSILON);
		ts_deboornet_free(&net);
		free(result);
		result = NULL;

		/* Tangent */
		C(ts_bspline_eval(&derivative, knots[i], &net, &status))
		C(ts_deboornet_result(&net, &result, &status))
		ts_vec_norm(result, 3, result);
		dist = ts_distance(result, frames[i].tangent, 3);
		CuAssertDblEquals(tc, 0, dist, POINT_EPSILON);
		ts_deboornet_free(&net);
		free(result);
		result = NULL;

		/* Normal */
		angle = (tsReal) 0.0;
		angle = ts_vec_angle(frames[i].normal,
		                     frames[i].tangent,
		                     NULL,
		                     3);
		CuAssertDblEquals(tc, 90.0, angle, ANGLE_EPSILON);
		angle = (tsReal) 0.0;
		angle = ts_vec_angle(frames[i].normal,
		                     frames[i].binormal,
		                     NULL,
		                     3);
		CuAssertDblEquals(tc, 90.0, angle, ANGLE_EPSILON);

		/* Binormal */
		angle = (tsReal) 0.0;
		angle = ts_vec_angle(frames[i].binormal,
		                     frames[i].tangent,
		                     NULL,
		                     3);
		CuAssertDblEquals(tc, 90.0, angle, ANGLE_EPSILON);
		angle = (tsReal) 0.0;
		angle = ts_vec_angle(frames[i].binormal,
		                     frames[i].normal,
		                     NULL,
		                     3);
		CuAssertDblEquals(tc, 90.0, angle, ANGLE_EPSILON);
	}

	___TEARDOWN___
	ts_bspline_free(&spline);
	ts_bspline_free(&derivative);
	ts_deboornet_free(&net);
	free(result);
}

#define TUBULAR_SEGMENTS 50
#define RADIAL_SEGMENTS 8

/* NOTE: If test_copies is true, some frames will repeat. */
void
rmf_test_tubular_mesh(CuTest *tc, tsReal vec_eps, tsReal dist_eps,
					  tsReal angle_eps, tsBSpline spline, size_t radius,
					  int test_copies) {
	___SETUP___
	tsReal knots[TUBULAR_SEGMENTS+1], v[3], w[3], line[3],
	       vertices[(TUBULAR_SEGMENTS + 1) * (RADIAL_SEGMENTS+1) * 3],
	       mag, t, normals[(TUBULAR_SEGMENTS+1) * (RADIAL_SEGMENTS+1) * 3],
		   tangents[(TUBULAR_SEGMENTS+1) * (RADIAL_SEGMENTS+1) * 4],
		   uvs[(TUBULAR_SEGMENTS+1) * (RADIAL_SEGMENTS+1) * 2];
	tsFrame frames[TUBULAR_SEGMENTS + 1 + 3];  /* + 3 in case of test_copies */
	size_t f_i, f_i_prev, i, j, k, a, b, c, d, x;
	size_t tubular_segments = TUBULAR_SEGMENTS;
	size_t indices[TUBULAR_SEGMENTS * RADIAL_SEGMENTS * 6];

	___GIVEN___
	ts_bspline_uniform_knot_seq(&spline, TUBULAR_SEGMENTS+1, knots);
	C(ts_bspline_compute_rmf(&spline,
							 knots,
							 TUBULAR_SEGMENTS+1,
							 0,
							 frames,
							 &status))
	if (test_copies) {
		/* We will repeat the first frame, the second frame, and the last frame.
		 * For this, we first have to shift the frames array to make space
		 * for the repeated frames. */
		memmove(frames + 4, frames + 2,
				sizeof(tsFrame) * (TUBULAR_SEGMENTS + 1 - 2));
		frames[3] = frames[2] = frames[1]; /* repeat second frame */
		frames[1] = frames[0]; /* repeat first frame */
		frames[TUBULAR_SEGMENTS + 1 + 2] = frames[TUBULAR_SEGMENTS + 1];
	}

	___WHEN___
	C(ts_frame_create_tube(frames, &tubular_segments, RADIAL_SEGMENTS, radius,
						   vertices, normals, tangents, uvs, indices, &status))

	___THEN___
	CuAssert(tc, "tubular_segments must not increase",
		     tubular_segments <= TUBULAR_SEGMENTS);
	for (i = 0; i <= tubular_segments; i++) {
		/* If some frames are repeated, we need to adjust
		 * the frames index accordingly. Note that this is not necessary
		 * in normal API usage, we just need to do so here because
		 * we need to compare the generated values with the
		 * values present in the frames array.*/
		if (!test_copies) {
			f_i = i;
		} else switch (i) {
			case 0: f_i = 0; break;
			case 1: f_i = 2; break;
			default: f_i = i + 2;
		}

		for (j = 0; j <= RADIAL_SEGMENTS; j++) {

			/* Check: Tangents must be the same. */
			for (k = 0; k < 3; k++) {
				CuAssertDblEquals(tc, frames[f_i].tangent[k],
								  tangents[i * (RADIAL_SEGMENTS+1) * 4
								           + j * 4 + k],
								  POINT_EPSILON);
			}
			CuAssertDblEquals(tc, 1.0, tangents[i * (RADIAL_SEGMENTS+1) * 4
			                                    + j * 4 + 3], POINT_EPSILON);

			/* Check: Normals must be unit vectors. */
			mag = ts_vec_mag(normals + i * (RADIAL_SEGMENTS+1) * 3 + j * 3, 3);
			CuAssertDblEquals(tc, 1.0, mag, POINT_EPSILON);

			if (i >= 1 && j >= 1) {
				/* Note that ABCD describes the rectangle counter-clockwise. */
				a = indices[(i - 1) * 6 * RADIAL_SEGMENTS
							+ (j - 1) * 6] * 3;
				b = indices[(i - 1) * 6 * RADIAL_SEGMENTS
							+ (j - 1) * 6 + 2] * 3;
				c = indices[(i - 1) * 6 * RADIAL_SEGMENTS
							+ (j - 1) * 6 + 5] * 3;
				d = indices[(i - 1) * 6 * RADIAL_SEGMENTS
							+ (j - 1) * 6 + 1] * 3;

				CuAssertIntEquals(tc, b, indices[(i - 1) * 6 * RADIAL_SEGMENTS
				       	                         + (j - 1) * 6 + 3] * 3);
				CuAssertIntEquals(tc, d, indices[(i - 1) * 6 * RADIAL_SEGMENTS
				                                 + (j - 1) * 6 + 4] * 3);

				if (i >= 2) {
					/* Check: This face should connect to previous face
					 * (this AD with last BC) */
					CuAssertIntEquals(tc, a,
									  indices[(i - 2) * 6 * RADIAL_SEGMENTS
									          + (j - 1) * 6 + 2] * 3);
					CuAssertIntEquals(tc, d,
									  indices[(i - 2) * 6 * RADIAL_SEGMENTS
									          + (j - 1) * 6 + 5] * 3);
				}

				/* The vector for the line drawn between frames. */
				ts_vec_sub(frames[f_i].position, frames[f_i_prev].position, 3,
						   line);

				/* Check: The three vectors AB, CD, and line must be equal. */
				ts_vec_sub(vertices + b, vertices + a, 3, v);
				CuAssertDblEquals(tc, 0.0, ts_distance(v, line, 3), vec_eps);
				ts_vec_sub(vertices + c, vertices + d, 3, v);
				CuAssertDblEquals(tc, 0.0, ts_distance(v, line, 3), vec_eps);

				/* Check: Distance from face corners (A,B,C,D) to line
				 * must equal radius. */
				for (k = 0; k < 4; k++) {
					switch (k) {
						case 0: x = a; break;
						case 1: x = b; break;
						case 2: x = d; break;
						case 3: x = c; break;
						default: CuFail(tc, "unreachable");
					}
					ts_vec_norm(line, 3, w);
					ts_vec_sub(frames[k % 2 ? f_i : f_i_prev].position,
							   vertices + x, 3, v);
					t = ts_vec_dot(v, w, 3);
					/* project corner onto line */
					ts_vec_mul(w, 3, t, w);
					ts_vec_add(frames[k % 2 ? f_i : f_i_prev].position,
							   w, 3, w);
					CuAssertDblEquals(tc, radius,
									  ts_distance(w, vertices + x, 3),
									  dist_eps);
				}

				ts_vec_add(vertices + a, vertices + c, 3, v);
				ts_vec_mul(v, 3, 0.5, v);
				ts_vec_add(vertices + b, vertices + d, 3, w);
				ts_vec_mul(w, 3, 0.5, w);
				/* Check: Center of AC and BD must be equal. */
				CuAssertDblEquals(tc, 0.0, ts_distance(v, w, 3), vec_eps);
				ts_vec_add(frames[f_i_prev].position, frames[f_i].position, 3,
						   w);
				ts_vec_mul(w, 3, 0.5, w);
				ts_vec_sub(v, w, 3, w);
				ts_vec_norm(w, 3, w);
				ts_vec_norm(line, 3, line);
				/* Check: 90-degree angle between line and face. */
				CuAssertDblEquals(tc, 90.0, ts_vec_angle(w, line, NULL, 3),
								  angle_eps);
			}

			/* TODO: Test UVs for correctness */
		}

		f_i_prev = f_i;
	}

	___TEARDOWN___
}

void
rmf_tubular_mesh_cubic_spline(CuTest *tc)
{
	___SETUP___
	size_t i, j;
	tsBSpline spline = ts_bspline_init();

	___GIVEN___
	C(ts_bspline_new_with_control_points(
			5, 3, 3, TS_CLAMPED, &spline, &status,
			100.0, 200.0, 0.0,   /* P1 */
			150.0, 220.0, 10.0,  /* P2 */
			190.0, 120.0, 50.0,  /* P3 */
			260.0, 70.0,  30.0,  /* P4 */
			300.0, 200.0, 20.0)) /* P5 */

	___WHEN___ ___THEN___
	for (i = 1; i <= 10; i++) {
		/* Tubes generated for non-linear splines may deviate from the tests
		 * in some ways, as the faces will be at times longer or shorter
		 * than the lines they are built around. Hence, the epsilons
		 * here are larger than for rmf_tubular_mesh_line.
		 * Additionally, a higher radius leads to decreased accuracy,
		 * thus we need to make the epsilons dependent on the radius. */
		rmf_test_tubular_mesh(tc, 0.4 + i * 0.2, 0.1 + i * 0.05, 2.5,
							  spline, (tsReal) i, i == 1);
	}

	___TEARDOWN___
	ts_bspline_free(&spline);
}

void
rmf_tubular_mesh_line(CuTest *tc)
{
	___SETUP___
	size_t i = 1;
	tsBSpline spline = ts_bspline_init();
	tsReal eps = ANGLE_EPSILON;
#ifdef TINYSPLINE_FLOAT_PRECISION
	/* ANGLE_EPSILON is too narrow with float precision,
	 * so we change the original margin of 1e-3f to 1e-2f. */
	eps = (tsReal) 1e-2f; /* from 1e-3f (0.001) */
#endif

	___GIVEN___
	C(ts_bspline_new_with_control_points(
			2, 3, 1, TS_CLAMPED, &spline, &status,
			100.0, 200.0, 0.0,   /* P1 */
			200.0, 240.0, 20.0)) /* P2 */

	___WHEN___ ___THEN___
	for (i = 1; i <= 10; i++) {
		rmf_test_tubular_mesh(tc, POINT_EPSILON, POINT_EPSILON, eps, spline,
							  (tsReal) i, i == 1);
	}

	___TEARDOWN___
	ts_bspline_free(&spline);
}

CuSuite *
get_rmf_suite()
{
	CuSuite* suite = CuSuiteNew();
	SUITE_ADD_TEST(suite, rmf_vectors_of_frames);
	SUITE_ADD_TEST(suite, rmf_tubular_mesh_cubic_spline);
	SUITE_ADD_TEST(suite, rmf_tubular_mesh_line);
	return suite;
}
