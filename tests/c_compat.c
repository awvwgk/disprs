#include <dftd3.h>
#include <dftd4.h>

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

_Static_assert(_Generic((dftd3_error)0, dftd4_error: 1, default: 0), "shared error type");
_Static_assert(_Generic((dftd3_structure)0, dftd4_structure: 1, default: 0),
               "shared structure type");

static void close_value(double actual, double expected, double tolerance) {
    if (!isfinite(actual) || fabs(actual - expected) > tolerance) {
        fprintf(stderr, "%.17g != %.17g (tolerance %.3g)\n", actual, expected, tolerance);
        assert(0);
    }
}

static void test_d3(void) {
    dftd3_error error = dftd3_new_error();
    const int numbers[] = {6, 6};
    const double positions[] = {0.0, 0.0, 0.0, 6.0, 0.0, 0.0};
    dftd3_structure structure = dftd3_new_structure(error, 2, numbers, positions, NULL, NULL);
    dftd3_model model = dftd3_new_d3_model(error, structure);
    dftd3_param param = dftd3_load_rational_damping(error, "pbe", false);
    double energy;
    dftd3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    assert(!dftd3_check_error(error));
    assert(fabs(energy - -0.0005341413931338267) < 1e-15);
    dftd3_delete(param);
    dftd3_delete(model);
    dftd3_delete(structure);
    dftd3_delete(error);
    assert(!param && !model && !structure && !error);
    dftd3_delete(error);
    dftd3_delete_error(NULL);
    dftd3_delete_structure(NULL);
    dftd3_delete_model(NULL);
    dftd3_delete_param(NULL);
    dftd3_delete_gcp(NULL);

    error = dftd3_new_error();
    assert(dftd3_get_version() == 10600);
    assert(dftd3_check_error(NULL));
    assert(!dftd3_has_feature("mpi") && !dftd3_has_feature("unknown"));
    const int atoms[] = {6, 8, 7};
    const double xyz[] = {0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0};
    structure = dftd3_new_structure(error, 3, atoms, xyz, NULL, NULL);
    model = dftd3_new_d3_model(error, structure);
    double cn[3], c6[9], cn_only[3], c6_only[9];
    disprs_d3_get_properties(error, structure, model, cn, c6);
    assert(!dftd3_check_error(error));
    disprs_d3_get_properties(error, structure, model, cn_only, NULL);
    disprs_d3_get_properties(error, structure, model, NULL, c6_only);
    disprs_d3_get_properties(error, structure, model, NULL, NULL);
    assert(!dftd3_check_error(error));
    for (int first = 0; first < 3; ++first) {
        close_value(cn[first], cn_only[first], 0.);
        assert(cn[first] > 0.);
        for (int second = 0; second < 3; ++second) {
            close_value(c6[3 * first + second], c6_only[3 * first + second], 0.);
            close_value(c6[3 * first + second], c6[3 * second + first], 0.);
            assert(c6[3 * first + second] > 0.);
        }
    }
    cn_only[0] = c6_only[0] = 42.;
    disprs_d3_get_properties(error, NULL, model, cn_only, c6_only);
    assert(dftd3_check_error(error) && cn_only[0] == 42. && c6_only[0] == 42.);
    disprs_d3_get_properties(error, structure, NULL, cn_only, c6_only);
    assert(dftd3_check_error(error) && cn_only[0] == 42. && c6_only[0] == 42.);
    disprs_d3_get_properties(error, structure, model, cn_only, c6_only);
    assert(!dftd3_check_error(error));
    dftd3_param explicit_params[] = {
        dftd3_new_zero_damping(error, 1., 1., 1., 1., 1., 14.),
        dftd3_new_rational_damping(error, 1., 1., 1., .4, 4., 14.),
        dftd3_new_mzero_damping(error, 1., 1., 1., 1., 1., 14., .1),
        dftd3_new_mrational_damping(error, 1., 1., 1., .4, 4., 14.),
        dftd3_new_optimizedpower_damping(error, 1., 1., 1., .4, 4., 14., 2.),
        dftd3_new_cso_damping(error, 1., 1., .86, 2.5, 0., 6.25, 14.),
        dftd3_new_z_damping(error, 1., 1., 1., 200770., 14.)};
    dftd3_param (*loaders[])(dftd3_error, char *, bool) = {
        dftd3_load_zero_damping,      dftd3_load_rational_damping,       dftd3_load_mzero_damping,
        dftd3_load_mrational_damping, dftd3_load_optimizedpower_damping, dftd3_load_cso_damping,
        dftd3_load_z_damping};
    for (int family = 0; family < 7; ++family) {
        double gradient[9], virial[9], pair2[9], pair3[9], hessian[81];
        param = explicit_params[family];
        assert(param && !dftd3_check_error(error));
        dftd3_get_dispersion(error, structure, model, param, &energy, gradient, virial);
        dftd3_get_pairwise_dispersion(error, structure, model, param, pair2, pair3);
        double sum = 0.;
        for (int index = 0; index < 9; ++index)
            sum += pair2[index] + pair3[index];
        close_value(sum, energy, 1e-13);
        dftd3_get_dispersion_hessian(error, structure, model, param, &sum, hessian);
        close_value(sum, energy, 1e-13);
        for (int column = 0; column < 9; ++column) {
            double displaced[9], plus[9], minus[9], sample;
            memcpy(displaced, xyz, sizeof xyz);
            displaced[column] += 1e-5;
            dftd3_update_structure(error, structure, displaced, NULL);
            dftd3_get_dispersion(error, structure, model, param, &sample, plus, NULL);
            displaced[column] -= 2e-5;
            dftd3_update_structure(error, structure, displaced, NULL);
            dftd3_get_dispersion(error, structure, model, param, &sample, minus, NULL);
            for (int row = 0; row < 9; ++row) {
                close_value(hessian[column * 9 + row], (plus[row] - minus[row]) / 2e-5, 1e-7);
                close_value(hessian[column * 9 + row], hessian[row * 9 + column], 1e-12);
            }
        }
        dftd3_update_structure(error, structure, xyz, NULL);
        assert(!dftd3_check_error(error));
        dftd3_delete(param);
        param = loaders[family](error, "pbe", true);
        if (family == 6) {
            assert(!param && dftd3_check_error(error));
        } else {
            assert(param && !dftd3_check_error(error));
            dftd3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
            assert(isfinite(energy) && !dftd3_check_error(error));
        }
        dftd3_delete(param);
    }

    dftd3_gcp gcp = dftd3_load_gcp_param(error, structure, "pbeh3c", NULL);
    assert(gcp && !dftd3_check_error(error));
    double eta, recovered_eta, slater[3], recovered_slater[3];
    disprs_d3_get_gcp_controls(error, gcp, &eta, NULL, NULL);
    disprs_d3_get_gcp_parameters(error, gcp, 3, NULL, NULL, NULL, NULL, NULL, slater, NULL, NULL);
    const double overflow_eta = DBL_MAX;
    disprs_d3_set_gcp_controls(error, gcp, &overflow_eta, NULL, NULL);
    assert(dftd3_check_error(error));
    disprs_d3_get_gcp_controls(error, gcp, &recovered_eta, NULL, NULL);
    assert(!dftd3_check_error(error) && eta == recovered_eta);
    disprs_d3_set_gcp_controls(error, gcp, &overflow_eta, NULL, NULL);
    assert(dftd3_check_error(error));
    disprs_d3_get_gcp_parameters(error, gcp, 3, NULL, NULL, NULL, NULL, NULL, recovered_slater,
                                 NULL, NULL);
    assert(!dftd3_check_error(error));
    for (int index = 0; index < 3; ++index) {
        assert(slater[index] == recovered_slater[index]);
    }
    disprs_d3_set_gcp_controls(error, gcp, &overflow_eta, NULL, NULL);
    assert(dftd3_check_error(error));
    disprs_d3_set_gcp_controls(error, gcp, &eta, NULL, NULL);
    assert(!dftd3_check_error(error));
    double gradient[9], virial[9], only_gradient[9], only_virial[9];
    dftd3_get_counterpoise(error, structure, gcp, &energy, gradient, virial);
    dftd3_get_counterpoise(error, structure, gcp, &energy, only_gradient, NULL);
    dftd3_get_counterpoise(error, structure, gcp, &energy, NULL, only_virial);
    for (int index = 0; index < 9; ++index) {
        close_value(gradient[index], only_gradient[index], 1e-15);
        close_value(virial[index], only_virial[index], 1e-15);
    }
    dftd3_set_gcp_mpi_comm(error, gcp, 0);
    assert(dftd3_check_error(error));
    dftd3_set_gcp_work_partition(error, gcp, 0, 1);
    assert(!dftd3_check_error(error));
    dftd3_set_model_mpi_comm(error, model, 0);
    assert(dftd3_check_error(error));
    dftd3_set_model_work_partition(error, model, 0, 1);
    assert(!dftd3_check_error(error));
    dftd3_set_model_realspace_cutoff(error, model, NAN, 20., 20.);
    assert(dftd3_check_error(error));
    dftd3_set_model_realspace_cutoff_smooth(error, model, 20., 20., 20., 1., 1.);
    assert(!dftd3_check_error(error));
    int invalid_ghost[] = {0, 3};
    dftd3_set_model_ghost_index(error, model, invalid_ghost, 2);
    assert(dftd3_check_error(error));
    param = dftd3_load_rational_damping(error, "pbe", false);
    dftd3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    double baseline = energy, bad_xyz[9];
    memcpy(bad_xyz, xyz, sizeof xyz);
    bad_xyz[2] = NAN;
    dftd3_update_structure(error, structure, bad_xyz, NULL);
    assert(dftd3_check_error(error));
    dftd3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    close_value(energy, baseline, 0.);
    assert(!dftd3_check_error(error));
    dftd3_param bad = dftd3_new_rational_damping(error, NAN, 1., 1., .4, 4., 14.);
    assert(!bad && dftd3_check_error(error));
    char message[512];
    memset(message, 'x', sizeof message);
    int capacity = 4;
    dftd3_get_error(error, message, &capacity);
    assert(message[3] == '\0' && message[4] == 'x');
    dftd3_get_error(error, message, NULL);
    assert(strlen(message) > 3);
    capacity = 0;
    message[0] = 'x';
    dftd3_get_error(error, message, &capacity);
    assert(message[0] == 'x');
    dftd3_structure smaller = dftd3_new_structure(error, 2, atoms, xyz, NULL, NULL);
    energy = 42.;
    dftd3_get_dispersion(error, smaller, model, param, &energy, NULL, NULL);
    assert(dftd3_check_error(error) && energy == 42.);
    dftd3_get_counterpoise(error, smaller, gcp, &energy, NULL, NULL);
    assert(dftd3_check_error(error) && energy == 42.);
    dftd3_delete(smaller);
    dftd3_delete(gcp);
    dftd3_delete(param);
    dftd3_delete(model);
    dftd3_delete(structure);
    dftd3_delete(error);
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters): independent axes of the test matrix.
static void test_d4(int kind, int dimensions) {
    dftd4_error error = dftd4_new_error();
    const int numbers[] = {6, 8, 1};
    const double positions[] = {.4, .8, 1.2, 4.4, 2.8, 3.6, 1.2, 4.5, 2.1};
    const double lattice[] = {12., 0., 0., .5, 13., 0., .3, .7, 14.};
    const bool periodic[] = {dimensions > 0, dimensions > 1, dimensions > 2};
    const double charge = .25;
    dftd4_structure mol = dftd4_new_structure(error, 3, numbers, positions, &charge,
                                              dimensions ? lattice : NULL, periodic);
    dftd4_model model;
    switch (kind) {
    case 0:
        model = dftd4_new_d4_model(error, mol);
        break;
    case 1:
        model = dftd4_new_d4s_model(error, mol);
        break;
    case 2:
        model = dftd4_custom_d4_model(error, mol, 2., 1., 4.);
        break;
    default:
        model = dftd4_custom_d4s_model(error, mol, 2., 1.);
        break;
    }
    assert(mol && model && !dftd4_check_error(error));
    dftd4_set_model_realspace_cutoff(error, model, 10., 8., 7.);
    dftd4_set_model_realspace_cutoff_smooth(error, model, 10., 8., 7., 1., 1.);
    dftd4_param param = dftd4_load_rational_damping(error, "pbe", true);
    double energy, scalar, gradient[9], virial[9], separate[9], pair2[9], pair3[9];
    double cn[3], charges[3], c6[9], alpha[3];
    dftd4_get_properties(error, mol, model, cn, charges, c6, alpha);
    close_value(charges[0] + charges[1] + charges[2], charge, 1e-12);
    for (int atom = 0; atom < 3; ++atom) {
        assert(isfinite(cn[atom]) && cn[atom] >= 0. && alpha[atom] > 0.);
        for (int other = 0; other < 3; ++other) {
            assert(c6[atom * 3 + other] > 0.);
            close_value(c6[atom * 3 + other], c6[other * 3 + atom], 1e-12);
        }
    }
    dftd4_get_properties(error, mol, model, NULL, NULL, separate, NULL);
    for (int index = 0; index < 9; ++index)
        close_value(separate[index], c6[index], 0.);
    dftd4_get_dispersion(error, mol, model, param, &energy, gradient, virial);
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, NULL);
    close_value(scalar, energy, 1e-14);
    const double original_scalar = scalar;
    dftd4_get_dispersion(error, mol, model, param, &scalar, separate, NULL);
    for (int index = 0; index < 9; ++index)
        close_value(separate[index], gradient[index], 1e-14);
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, separate);
    for (int index = 0; index < 9; ++index)
        close_value(separate[index], virial[index], 1e-14);
    dftd4_get_pairwise_dispersion(error, mol, model, param, pair2, pair3);
    double sum = 0.;
    for (int index = 0; index < 9; ++index)
        sum += pair2[index] + pair3[index];
    close_value(sum, energy, 1e-13);
    dftd4_param explicit_param =
        dftd4_new_rational_damping(error, 1., .95948085, 1., .38574991, 4.80688534, 16.);
    dftd4_get_dispersion(error, mol, model, explicit_param, &scalar, NULL, NULL);
    close_value(scalar, energy, 1e-14);
    dftd4_delete(explicit_param);
    double storage[83];
    storage[0] = storage[82] = 12345.;
    double *hessian = storage + 1;
    dftd4_get_numerical_hessian(error, mol, model, param, hessian);
    assert(storage[0] == 12345. && storage[82] == 12345.);
    assert(!dftd4_check_error(error));
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, NULL);
    close_value(scalar, original_scalar, 0.);
    for (int column = 0; column < 9; ++column) {
        double displaced[9], plus[9], minus[9], eplus, eminus;
        memcpy(displaced, positions, sizeof positions);
        displaced[column] += 2e-4;
        dftd4_update_structure(error, mol, displaced, NULL);
        dftd4_get_dispersion(error, mol, model, param, &eplus, plus, NULL);
        displaced[column] -= 4e-4;
        dftd4_update_structure(error, mol, displaced, NULL);
        dftd4_get_dispersion(error, mol, model, param, &eminus, minus, NULL);
        close_value(gradient[column], (eplus - eminus) / 4e-4, 2e-8);
        for (int row = 0; row < 9; ++row) {
            close_value(hessian[column * 9 + row], (plus[row] - minus[row]) / 4e-4, 2e-7);
            close_value(hessian[column * 9 + row], hessian[row * 9 + column], 2e-7);
        }
        double samples[2];
        for (int side = 0; side < 2; ++side) {
            double strained[9];
            const double step = (side ? 1. : -1.) * 1e-5;
            memcpy(displaced, positions, sizeof positions);
            memcpy(strained, lattice, sizeof lattice);
            for (int atom = 0; atom < 3; ++atom) {
                displaced[3 * atom + column % 3] += step * positions[3 * atom + column / 3];
                strained[3 * atom + column % 3] += step * lattice[3 * atom + column / 3];
            }
            dftd4_update_structure(error, mol, displaced, dimensions ? strained : NULL);
            dftd4_get_dispersion(error, mol, model, param, &samples[side], NULL, NULL);
        }
        close_value(virial[column], (samples[1] - samples[0]) / 2e-5, 2e-8);
        dftd4_update_structure(error, mol, positions, dimensions ? lattice : NULL);
    }
    assert(!dftd4_check_error(error));
    double bad[9];
    memcpy(bad, positions, sizeof positions);
    bad[0] = NAN;
    dftd4_update_structure(error, mol, bad, NULL);
    assert(dftd4_check_error(error));
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, NULL);
    close_value(scalar, original_scalar, 0.);
    assert(!dftd4_check_error(error));
    dftd4_set_model_realspace_cutoff(error, model, NAN, 8., 7.);
    assert(dftd4_check_error(error));
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, NULL);
    close_value(scalar, original_scalar, 0.);
    dftd4_set_model_realspace_cutoff(error, model, 1., 1., 7.);
    dftd4_get_dispersion(error, mol, model, param, &scalar, NULL, NULL);
    assert(fabs(scalar - energy) > 1e-8);
    dftd4_delete(param);
    dftd4_delete(model);
    dftd4_delete(mol);
    dftd4_delete(error);
    assert(!param && !model && !mol && !error);
}

static void test_d4_errors(void) {
    assert(dftd4_get_version() == 40200 && dftd4_check_error(NULL));
    dftd4_error error = dftd4_new_error();
    double output = 42.;
    dftd4_get_dispersion(error, NULL, NULL, NULL, &output, NULL, NULL);
    assert(dftd4_check_error(error) && output == 42.);
    dftd4_get_properties(error, NULL, NULL, &output, NULL, NULL, NULL);
    assert(dftd4_check_error(error) && output == 42.);
    dftd4_get_pairwise_dispersion(error, NULL, NULL, NULL, &output, &output);
    assert(dftd4_check_error(error) && output == 42.);
    dftd4_get_numerical_hessian(error, NULL, NULL, NULL, &output);
    assert(dftd4_check_error(error) && output == 42.);
    assert(!dftd4_new_d4_model(error, NULL) && dftd4_check_error(error));
    dftd4_update_structure(error, NULL, &output, NULL);
    assert(dftd4_check_error(error));
    dftd4_set_model_realspace_cutoff(error, NULL, 10., 8., 7.);
    assert(dftd4_check_error(error));
    assert(!dftd4_load_rational_damping(error, "unknown-method", true));
    assert(dftd4_check_error(error));
    assert(!dftd4_new_rational_damping(error, NAN, 1., 1., .4, 4., 16.));
    assert(dftd4_check_error(error));
    char message[512];
    dftd4_get_error(error, message, NULL);
    assert(strlen(message) > 0);
    const int numbers[] = {6, 8};
    const double xyz[] = {0., 0., 0., 4., 1., 2.};
    const double infinite_charge = INFINITY;
    assert(!dftd4_new_structure(error, 2, numbers, xyz, &infinite_charge, NULL, NULL));
    assert(dftd4_check_error(error));
    dftd4_structure mol = dftd4_new_structure(error, 2, numbers, xyz, NULL, NULL, NULL);
    assert(mol && !dftd4_check_error(error));
    assert(!dftd4_custom_d4_model(error, mol, NAN, 2., 6.));
    assert(dftd4_check_error(error));
    dftd4_model model = dftd4_new_d4_model(error, mol);
    dftd4_param param = dftd4_load_rational_damping(error, "pbe", false);
    dftd4_get_numerical_hessian(error, mol, model, param, NULL);
    assert(dftd4_check_error(error));
    dftd4_get_dispersion(error, mol, model, param, &output, NULL, NULL);
    assert(!dftd4_check_error(error));
    dftd4_delete(param);
    dftd4_delete(model);
    dftd4_delete(mol);
    dftd4_delete(error);
    dftd4_delete(error);
    dftd4_delete_error(NULL);
    dftd4_delete_structure(NULL);
    dftd4_delete_model(NULL);
    dftd4_delete_param(NULL);
}

static void test_periodic_defaults(void) {
    const int numbers[] = {6, 8};
    const double xyz[] = {0., 0., 0., 4., 1., 2.};
    const double lattice[] = {8., 0., 0., 0., 9., 0., 0., 0., 10.};
    const bool periodic[] = {true, true, true};
    const bool molecular[] = {false, false, false};
    double d3_energy[3], d4_energy[3];
    dftd3_error error3 = dftd3_new_error();
    dftd4_error error4 = dftd4_new_error();
    for (int mode = 0; mode < 3; ++mode) {
        const bool *flags = mode == 0 ? NULL : (mode == 1 ? periodic : molecular);
        dftd3_structure mol3 = dftd3_new_structure(error3, 2, numbers, xyz, lattice, flags);
        dftd3_model model3 = dftd3_new_d3_model(error3, mol3);
        dftd3_param param3 = dftd3_load_rational_damping(error3, "pbe", false);
        dftd3_set_model_realspace_cutoff(error3, model3, 12., 10., 10.);
        dftd3_get_dispersion(error3, mol3, model3, param3, &d3_energy[mode], NULL, NULL);
        assert(!dftd3_check_error(error3));
        dftd4_structure mol4 = dftd4_new_structure(error4, 2, numbers, xyz, NULL, lattice, flags);
        dftd4_model model4 = dftd4_new_d4_model(error4, mol4);
        dftd4_param param4 = dftd4_load_rational_damping(error4, "pbe", false);
        dftd4_set_model_realspace_cutoff(error4, model4, 12., 10., 10.);
        dftd4_get_dispersion(error4, mol4, model4, param4, &d4_energy[mode], NULL, NULL);
        assert(!dftd4_check_error(error4));
        dftd3_delete(param3);
        dftd3_delete(model3);
        dftd3_delete(mol3);
        dftd4_delete(param4);
        dftd4_delete(model4);
        dftd4_delete(mol4);
    }
    close_value(d3_energy[0], d3_energy[1], 0.);
    close_value(d4_energy[0], d4_energy[1], 0.);
    assert(fabs(d3_energy[0] - d3_energy[2]) > 1e-6);
    assert(fabs(d4_energy[0] - d4_energy[2]) > 1e-6);
    dftd3_delete(error3);
    dftd4_delete(error4);
}

static void test_nonfinite_results(void) {
    const int numbers[] = {6, 6};
    const double positions[] = {0., 0., 0., .1, 0., 0.};
    for (int kind = 3; kind <= 4; ++kind) {
        void *error = disprs_d3_new_error();
        void *structure =
            kind == 3 ? disprs_d3_new_structure(error, 2, numbers, positions, NULL, NULL)
                      : disprs_d4_new_structure(error, 2, numbers, positions, NULL, NULL, NULL);
        void *model = kind == 3 ? disprs_d3_new_model(error, structure)
                                : disprs_d4_new_model(error, structure, 0);
        void *param =
            kind == 3 ? disprs_d3_new_rational_damping(error, DBL_MAX, DBL_MAX, 0., 0., 0., 14.)
                      : disprs_d4_new_rational_damping(error, DBL_MAX, DBL_MAX, 0., 0., 0., 16.);
        assert(structure && model && param && !disprs_d3_check_error(error));
        double energy = 42., gradient[6], virial[9], pair2[4], pair3[4], hessian[36];
        for (int index = 0; index < 36; ++index)
            hessian[index] = 42.;
        for (int index = 0; index < 9; ++index)
            virial[index] = 42.;
        for (int index = 0; index < 6; ++index)
            gradient[index] = 42.;
        for (int index = 0; index < 4; ++index)
            pair2[index] = pair3[index] = 42.;
        void (*dispersion)(void *, void *, void *, void *, double *, double *, double *) =
            kind == 3 ? disprs_d3_get_dispersion : disprs_d4_get_dispersion;
        dispersion(error, structure, model, param, &energy, NULL, NULL);
        assert(disprs_d3_check_error(error) && energy == 42.);
        dispersion(error, structure, model, param, &energy, gradient, virial);
        assert(disprs_d3_check_error(error) && energy == 42.);
        if (kind == 3) {
            disprs_d3_get_pairwise_dispersion(error, structure, model, param, pair2, pair3);
            assert(disprs_d3_check_error(error));
            disprs_d3_get_dispersion_hessian(error, structure, model, param, &energy, hessian);
        } else {
            disprs_d4_get_pairwise_dispersion(error, structure, model, param, pair2, pair3);
            assert(disprs_d3_check_error(error));
            disprs_d4_get_numerical_hessian(error, structure, model, param, hessian);
        }
        assert(disprs_d3_check_error(error) && energy == 42.);
        for (int index = 0; index < 36; ++index)
            assert(hessian[index] == 42.);
        for (int index = 0; index < 9; ++index)
            assert(virial[index] == 42.);
        for (int index = 0; index < 6; ++index)
            assert(gradient[index] == 42.);
        for (int index = 0; index < 4; ++index)
            assert(pair2[index] == 42. && pair3[index] == 42.);
        if (kind == 3) {
            disprs_d3_delete_param(&param);
            param = disprs_d3_load_param(error, 1, "pbe", false);
        } else {
            disprs_d4_delete_param(&param);
            param = disprs_d4_load_param(error, "pbe", false);
        }
        dispersion(error, structure, model, param, &energy, gradient, virial);
        assert(!disprs_d3_check_error(error) && isfinite(energy));
        if (kind == 3) {
            disprs_d3_delete_param(&param);
            disprs_d3_delete_model(&model);
            disprs_d3_delete_structure(&structure);
        } else {
            disprs_d4_delete_param(&param);
            disprs_d4_delete_model(&model);
            disprs_d4_delete_structure(&structure);
        }
        disprs_d3_delete_error(&error);
    }
}

int main(void) {
    test_d3();
    test_d4_errors();
    test_nonfinite_results();
    test_periodic_defaults();
    for (int kind = 0; kind < 4; ++kind)
        for (int dimensions = 0; dimensions < 4; ++dimensions)
            test_d4(kind, dimensions);
    puts("D3/D4 C compatibility checks passed");
    return 0;
}