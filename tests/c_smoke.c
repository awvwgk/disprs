#include <disprs.h>

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

static void test_shared_handles(void) {
    const int numbers[] = {6, 6};
    const double positions[] = {0.0, 0.0, 0.0, 6.0, 0.0, 0.0};
    const double displaced[] = {0.0, 0.0, 0.0, 7.0, 0.0, 0.0};
    const double charge = 0.5;
    for (int origin = 0; origin < 3; ++origin) {
        disprs_error error = origin == 0   ? disprs_new_error()
                             : origin == 1 ? disprs_d3_new_error()
                                           : disprs_d4_new_error();
        disprs_structure structure =
            origin == 0 ? disprs_new_structure(error, 2, numbers, positions, &charge, NULL, NULL)
            : origin == 1
                ? disprs_d3_new_structure(error, 2, numbers, positions, NULL, NULL)
                : disprs_d4_new_structure(error, 2, numbers, positions, &charge, NULL, NULL);
        assert(structure && !disprs_check_error(error));
        disprs_d3_model d3 = disprs_d3_new_model(error, structure);
        disprs_d4_model d4 = disprs_d4_new_model(error, structure, DISPRS_D4);
        disprs_d3_param param3 = disprs_d3_load_param(error, DISPRS_D3_RATIONAL, "pbe", false);
        disprs_d4_param param4 = disprs_d4_load_param(error, "pbe", false);
        assert(d3 && d4 && param3 && param4 && !disprs_check_error(error));
        double energy, charges[2];
        disprs_d3_get_dispersion(error, structure, d3, param3, &energy, NULL, NULL);
        assert(!disprs_check_error(error) && fabs(energy + 0.0005341413931338267) < 1e-15);
        if (origin == 0)
            disprs_d4_update_structure(error, structure, displaced, NULL);
        if (origin == 1)
            disprs_update_structure(error, structure, displaced, NULL);
        if (origin == 2)
            disprs_d3_update_structure(error, structure, displaced, NULL);
        assert(!disprs_check_error(error));
        disprs_d4_get_properties(error, structure, d4, NULL, charges, NULL, NULL);
        assert(!disprs_check_error(error));
        assert(fabs(charges[0] + charges[1] - (origin == 1 ? 0.0 : charge)) < 1e-12);
        disprs_d4_get_dispersion(error, structure, d4, param4, &energy, NULL, NULL);
        assert(!disprs_check_error(error) && isfinite(energy));
        const double saved_energy = energy;
        disprs_update_structure(error, structure, NULL, NULL);
        assert(disprs_d3_check_error(error) && disprs_d4_check_error(error));
        char message[128];
        const int capacity = sizeof message;
        disprs_get_error(error, message, &capacity);
        assert(message[0]);
        disprs_d4_get_dispersion(error, structure, d4, param4, &energy, NULL, NULL);
        assert(!disprs_check_error(error) && energy == saved_energy);
        disprs_d3_gcp gcp = disprs_d3_load_gcp(error, structure, "pbeh3c", NULL);
        assert(gcp && !disprs_check_error(error));
        disprs_d3_get_counterpoise(error, structure, gcp, &energy, NULL, NULL);
        assert(!disprs_check_error(error) && isfinite(energy));
        disprs_d3_delete_gcp(&gcp);

        const int unsupported[] = {104, 104};
        disprs_structure other =
            disprs_new_structure(error, 2, unsupported, positions, NULL, NULL, NULL);
        assert(other && !disprs_check_error(error));
        assert(!disprs_d3_new_model(error, other) && disprs_check_error(error));
        assert(!disprs_d4_new_model(error, other, DISPRS_D4) && disprs_check_error(error));
        energy = 42.0;
        disprs_d3_get_dispersion(error, other, d3, param3, &energy, NULL, NULL);
        assert(disprs_check_error(error) && energy == 42.0);
        disprs_d4_get_dispersion(error, other, d4, param4, &energy, NULL, NULL);
        assert(disprs_check_error(error) && energy == 42.0);
        disprs_d3_update_structure(error, other, displaced, NULL);
        assert(!disprs_check_error(error));
        disprs_delete_structure(&other);
        disprs_d3_delete_param(&param3);
        disprs_d4_delete_param(&param4);
        disprs_d3_delete_model(&d3);
        disprs_d4_delete_model(&d4);
        if (origin == 0) {
            disprs_d3_delete_structure(&structure);
            disprs_d4_delete_error(&error);
        }
        if (origin == 1) {
            disprs_d4_delete_structure(&structure);
            disprs_delete_error(&error);
        }
        if (origin == 2) {
            disprs_delete_structure(&structure);
            disprs_d3_delete_error(&error);
        }
        assert(!structure && !error);
    }
    disprs_delete_structure(NULL);
    disprs_delete_error(NULL);
}

int main(void) {
    test_shared_handles();
    assert(strcmp(disprs_get_version(), "0.1.0") == 0);
    assert(disprs_d3_get_version() >= 10600);
    assert(!disprs_d3_has_feature("not-a-feature"));

    disprs_d3_error error = disprs_d3_new_error();
    const int numbers[] = {6, 6};
    const double positions[] = {0.0, 0.0, 0.0, 6.0, 0.0, 0.0};
    disprs_d3_structure structure =
        disprs_d3_new_structure(error, 2, numbers, positions, NULL, NULL);
    disprs_d3_model model = disprs_d3_new_model(error, structure);
    for (int kind = DISPRS_D3; kind <= DISPRS_D3S; ++kind) {
        disprs_d3_model selected = disprs_d3_new_model_kind(error, structure, kind);
        assert(selected && !disprs_d3_check_error(error));
        double selected_cn[2], selected_c6[4], legacy_cn[2], legacy_c6[4];
        disprs_d3_model legacy = kind == DISPRS_D3 ? disprs_d3_new_model(error, structure)
                                                   : disprs_d3_new_smooth_model(error, structure);
        disprs_d3_get_properties(error, structure, selected, selected_cn, selected_c6);
        assert(!disprs_d3_check_error(error));
        disprs_d3_get_properties(error, structure, legacy, legacy_cn, legacy_c6);
        assert(!disprs_d3_check_error(error));
        for (int index = 0; index < 2; ++index) {
            assert(selected_cn[index] == legacy_cn[index]);
        }
        for (int index = 0; index < 4; ++index) {
            assert(selected_c6[index] == legacy_c6[index]);
        }
        disprs_d3_delete_model(&selected);
        disprs_d3_delete_model(&legacy);
    }
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange): test rejection of an invalid
    // selector.
    assert(!disprs_d3_new_model_kind(error, structure, 2));
    assert(disprs_d3_check_error(error));
    assert(!disprs_d3_new_model_kind(error, NULL, DISPRS_D3S));
    assert(disprs_d3_check_error(error));
    disprs_d3_set_model_realspace_cutoff(error, model, 60.0, 40.0, 30.0);
    disprs_d3_set_model_realspace_cutoff_smooth(error, model, 60.0, 40.0, 30.0, 2.0, 2.0);
    disprs_d3_set_model_work_partition(error, model, 0, 1);
    disprs_d3_param param = disprs_d3_load_param(error, DISPRS_D3_RATIONAL, "pbe", false);
    double energy;
    double gradient[6];
    double virial[9];
    double pair2[4];
    double pair3[4];
    double hessian[36];
    double response_cn[2], response_c6[4], response_cartesian[12], response_strain[18];
    disprs_d3_get_property_response(error, structure, model, response_cn, response_c6,
                                    response_cartesian, response_strain, NULL, NULL);
    assert(!disprs_d3_check_error(error));
    assert(fabs(response_cartesian[0] + response_cartesian[3]) < 1e-12);
    assert(fabs(response_strain[0] - 6.0 * response_cartesian[3]) < 1e-12);
    response_cn[0] = response_cartesian[0] = response_strain[0] = 42.0;
    disprs_d3_get_property_response(error, structure, NULL, response_cn, NULL, response_cartesian,
                                    response_strain, NULL, NULL);
    assert(disprs_d3_check_error(error));
    assert(response_cn[0] == 42.0 && response_cartesian[0] == 42.0 && response_strain[0] == 42.0);
    disprs_d3_get_property_response(error, structure, model, NULL, NULL, NULL, NULL, NULL, NULL);
    assert(!disprs_d3_check_error(error));
    disprs_d3_get_dispersion(error, structure, model, param, &energy, gradient, virial);
    disprs_d3_get_pairwise_dispersion(error, structure, model, param, pair2, pair3);
    disprs_d3_get_dispersion_hessian(error, structure, model, param, &energy, hessian);

    assert(!disprs_d3_check_error(error));
    assert(isfinite(energy) && energy < 0.0);
    assert(fabs(gradient[0] + gradient[3]) < 1e-12);
    assert(isfinite(pair2[1]));
    assert(isfinite(hessian[0]));
    disprs_d3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    assert(fabs(energy - -0.0005341413931338267) < 1e-15);

    disprs_d3_param explicit_param =
        disprs_d3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
    assert(explicit_param != NULL);
    disprs_d3_get_dispersion(error, structure, model, explicit_param, &energy, gradient, virial);
    assert(!disprs_d3_check_error(error));
    assert(fabs(energy - -0.0005341413931338267) < 1e-15);
    assert(fabs(gradient[0] - -0.00018488441198083488) < 1e-14);
    assert(fabs(virial[0] - 0.0011093064718850092) < 1e-13);
    disprs_d3_get_dispersion_hessian(error, structure, model, explicit_param, &energy, hessian);
    assert(fabs(hessian[0] - 4.9042518143140056e-05) < 1e-13);
    assert(fabs(hessian[3] - -4.9042518143140056e-05) < 1e-13);
    disprs_d3_get_pairwise_dispersion(error, structure, model, explicit_param, pair2, pair3);
    assert(fabs(pair2[1] - -0.00026707069656691336) < 1e-15);
    assert(pair3[0] == 0.0 && pair3[1] == 0.0);
    const double smooth_positions[] = {0.0, 0.0, 0.0, 5.5, 0.0, 0.0};
    disprs_d3_update_structure(error, structure, smooth_positions, NULL);
    disprs_d3_set_model_realspace_cutoff_smooth(error, model, 6.0, 40.0, 30.0, 2.0, 0.0);
    disprs_d3_get_dispersion(error, structure, model, explicit_param, &energy, gradient, virial);
    assert(fabs(energy - -6.405046062946415e-5) < 1e-15);
    assert(fabs(gradient[0] - -3.4189782676699354e-4) < 1e-14);
    disprs_d3_update_structure(error, structure, positions, NULL);
    disprs_d3_set_model_realspace_cutoff(error, model, 60.0, 40.0, 30.0);
    disprs_d3_delete_param(&explicit_param);

    disprs_d3_gcp gcp = disprs_d3_load_gcp(error, structure, "pbeh3c", NULL);
    disprs_d3_set_gcp_realspace_cutoff(error, gcp, 30.0, 30.0);
    disprs_d3_set_gcp_work_partition(error, gcp, 0, 1);
    disprs_d3_get_counterpoise(error, structure, gcp, &energy, gradient, virial);
    disprs_d3_get_counterpoise_hessian(error, structure, gcp, &energy, hessian);
    assert(!disprs_d3_check_error(error));
    assert(fabs(energy - 2.581420786404301e-05) < 1e-15);
    assert(fabs(gradient[0] - 5.324722896549186e-05) < 1e-14);
    assert(fabs(virial[0] - -0.00031948337379295114) < 1e-13);
    assert(fabs(hessian[0] - 9.51416946673085e-05) < 1e-13);
    disprs_d3_delete_gcp(&gcp);

    disprs_d3_delete_param(&param);
    disprs_d3_delete_model(&model);
    disprs_d3_delete_structure(&structure);
    disprs_d3_delete_error(&error);

    error = disprs_d3_new_error();
    const int periodic_numbers[] = {6, 8};
    const double periodic_positions[] = {0.4, 0.8, 1.2, 4.4, 2.8, 3.6};
    const double lattice[] = {8.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 8.0};
    const bool periodic[] = {true, true, true};
    structure =
        disprs_d3_new_structure(error, 2, periodic_numbers, periodic_positions, lattice, periodic);
    model = disprs_d3_new_model(error, structure);
    disprs_d3_set_model_ewald(error, model, 10000, 1e-4, 10.0, -1);
    param = disprs_d3_new_rational_damping(error, 1.0, 0.7875, 0.0, 0.4289, 4.4407, 14.0);
    disprs_d3_get_dispersion(error, structure, model, param, &energy, gradient, virial);
    assert(!disprs_d3_check_error(error));
    assert(fabs(energy - -0.0028149228007689916) < 1e-12);
    assert(fabs(gradient[1] - 2.939668061386282e-5) < 1e-12);
    assert(fabs(gradient[2] - 5.2457766011066856e-5) < 1e-12);
    assert(fabs(virial[0] - 0.003746104731458539) < 1e-12);
    double partition_energy[2];
    disprs_d3_set_model_work_partition(error, model, 0, 2);
    disprs_d3_get_dispersion(error, structure, model, param, &partition_energy[0], NULL, NULL);
    disprs_d3_set_model_work_partition(error, model, 1, 2);
    disprs_d3_get_dispersion(error, structure, model, param, &partition_energy[1], NULL, NULL);
    assert(fabs(partition_energy[0] + partition_energy[1] - energy) < 1e-12);
    disprs_d3_set_model_work_partition(error, model, 0, 1);
    disprs_d3_set_model_ewald(error, model, 10000, 1e-4, 10.0, 64);
    disprs_d3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    assert(fabs(energy - -0.0028149228007021015) < 1e-11);
    disprs_d3_set_model_ewald(error, model, 10000, 1e-4, 10.0, 0);
    disprs_d3_get_dispersion(error, structure, model, param, &energy, NULL, NULL);
    assert(fabs(energy - -0.0028149227957823923) < 1e-10);
    gcp = disprs_d3_load_gcp(error, structure, "pbeh3c", NULL);
    disprs_d3_set_gcp_realspace_cutoff(error, gcp, 12.0, 12.0);
    disprs_d3_get_counterpoise(error, structure, gcp, &energy, gradient, virial);
    assert(fabs(energy - 0.0007503293440007774) < 1e-13);
    assert(fabs(gradient[1] - 0.0004289849315224243) < 1e-13);
    assert(fabs(gradient[2] - 0.0004974479116136508) < 1e-13);
    assert(fabs(virial[0] + 0.003558225718407879) < 1e-12);
    disprs_d3_get_counterpoise_hessian(error, structure, gcp, &energy, hessian);
    assert(fabs(energy - 0.0007503293440007774) < 1e-13);
    assert(fabs(hessian[0] - 0.0006505276577134137) < 1e-12);
    assert(fabs(hessian[3] + 0.0006505276577134137) < 1e-12);
    disprs_d3_delete_gcp(&gcp);
    disprs_d3_delete_param(&param);
    disprs_d3_delete_model(&model);
    disprs_d3_delete_structure(&structure);
    disprs_d3_delete_error(&error);

    error = disprs_d3_new_error();
    const int atm_numbers[] = {6, 8, 7};
    const double atm_positions[] = {0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 4.0, 0.0};
    structure = disprs_d3_new_structure(error, 3, atm_numbers, atm_positions, NULL, NULL);
    model = disprs_d3_new_model(error, structure);
    param = disprs_d3_new_rational_damping(error, 1.0, 0.7875, 1.0, 0.4289, 4.4407, 14.0);
    double atm_pair2[9];
    double atm_pair3[9];
    double atm_gradient[9];
    double atm_virial[9];
    disprs_d3_get_dispersion(error, structure, model, param, &energy, atm_gradient, atm_virial);
    disprs_d3_get_pairwise_dispersion(error, structure, model, param, atm_pair2, atm_pair3);
    double resolved = 0.0;
    for (int index = 0; index < 9; index++)
        resolved += atm_pair2[index] + atm_pair3[index];
    assert(!disprs_d3_check_error(error));
    assert(fabs(atm_pair3[1] - 2.962135990925684e-8) < 1e-17);
    assert(fabs(resolved - energy) < 1e-14);
    disprs_d3_delete_param(&param);
    disprs_d3_delete_model(&model);
    disprs_d3_delete_structure(&structure);
    disprs_d3_delete_error(&error);

    disprs_d4_error d4_error = disprs_d4_new_error();
    disprs_d4_structure d4_structure =
        disprs_d4_new_structure(d4_error, 2, numbers, positions, NULL, NULL, NULL);
    disprs_d4_model d4_model = disprs_d4_new_model(d4_error, d4_structure, DISPRS_D4);
    double fixed_charges[2] = {0.2, -0.3}, actual_charges[2], dqdr[12], dqdstrain[18];
    disprs_d4_set_fixed_charges(d4_error, d4_model, fixed_charges, 2);
    assert(!disprs_d4_check_error(d4_error));
    fixed_charges[0] = NAN;
    disprs_d4_set_fixed_charges(d4_error, d4_model, fixed_charges, 2);
    assert(disprs_d4_check_error(d4_error));
    disprs_d4_set_fixed_charges(d4_error, d4_model, fixed_charges, 1);
    assert(disprs_d4_check_error(d4_error));
    disprs_d4_get_property_response(d4_error, d4_structure, d4_model, NULL, actual_charges, NULL,
                                    NULL, NULL, NULL, dqdr, dqdstrain, NULL, NULL, NULL, NULL);
    assert(!disprs_d4_check_error(d4_error));
    assert(actual_charges[0] == 0.2 && actual_charges[1] == -0.3);
    for (int index = 0; index < 12; index++)
        assert(dqdr[index] == 0.0);
    for (int index = 0; index < 18; index++)
        assert(dqdstrain[index] == 0.0);
    disprs_d4_set_fixed_charges(d4_error, d4_model, NULL, 0);
    assert(!disprs_d4_check_error(d4_error));
    disprs_d4_get_properties(d4_error, d4_structure, d4_model, NULL, actual_charges, NULL, NULL);
    assert(fabs(actual_charges[0]) < 1e-12 && fabs(actual_charges[1]) < 1e-12);
    disprs_d4_get_property_response(d4_error, d4_structure, d4_model, response_cn, NULL,
                                    response_c6, NULL, response_cartesian, response_strain, NULL,
                                    NULL, NULL, NULL, NULL, NULL);
    assert(!disprs_d4_check_error(d4_error));
    assert(fabs(response_cartesian[0] + response_cartesian[3]) < 1e-12);
    assert(fabs(response_strain[0] - 6.0 * response_cartesian[3]) < 1e-12);
    response_cn[0] = response_cartesian[0] = response_strain[0] = 42.0;
    disprs_d4_get_property_response(d4_error, d4_structure, NULL, response_cn, NULL, NULL, NULL,
                                    response_cartesian, response_strain, NULL, NULL, NULL, NULL,
                                    NULL, NULL);
    assert(disprs_d4_check_error(d4_error));
    assert(response_cn[0] == 42.0 && response_cartesian[0] == 42.0 && response_strain[0] == 42.0);
    disprs_d4_get_property_response(d4_error, d4_structure, d4_model, NULL, NULL, NULL, NULL, NULL,
                                    NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    assert(!disprs_d4_check_error(d4_error));
    disprs_d4_param d4_param = disprs_d4_load_param(d4_error, "pbe", false);
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);

    assert(!disprs_d4_check_error(d4_error));
    assert(isfinite(energy) && energy < 0.0);
    assert(fabs(gradient[0] + gradient[3]) < 1e-12);
    double analytical_hessian[36], numerical_hessian[36];
    disprs_d4_get_dispersion_hessian(d4_error, d4_structure, d4_model, d4_param,
                                     analytical_hessian);
    assert(!disprs_d4_check_error(d4_error));
    disprs_d4_get_numerical_hessian(d4_error, d4_structure, d4_model, d4_param, numerical_hessian);
    assert(!disprs_d4_check_error(d4_error));
    for (int index = 0; index < 36; index++) {
        assert(fabs(analytical_hessian[index] - numerical_hessian[index]) < 1e-10);
        assert(fabs(analytical_hessian[index] - analytical_hessian[(index % 6) * 6 + index / 6]) <
               1e-12);
    }
    analytical_hessian[0] = 42.0;
    disprs_d4_get_dispersion_hessian(d4_error, d4_structure, d4_model, NULL, analytical_hessian);
    assert(disprs_d4_check_error(d4_error) && analytical_hessian[0] == 42.0);
    disprs_d4_get_dispersion_hessian(d4_error, d4_structure, d4_model, d4_param, NULL);
    assert(disprs_d4_check_error(d4_error));
    disprs_d4_get_dispersion_hessian(d4_error, d4_structure, d4_model, d4_param,
                                     analytical_hessian);
    assert(!disprs_d4_check_error(d4_error));
    double previous_energy = energy;
    (void)previous_energy;
    disprs_d4_delete_param(&d4_param);
    d4_param = disprs_d4_new_rational_damping(d4_error, 1.0, 0.95948085, 0.0, 0.38574991,
                                              4.80688534, 16.0);
    assert(d4_param != NULL);
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);
    assert(energy == previous_energy);
    const double updated_positions[] = {0.0, 0.0, 0.0, 6.1, 0.0, 0.0};
    disprs_d4_update_structure(d4_error, d4_structure, updated_positions, NULL);
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);
    assert(energy != previous_energy);

    disprs_d4_delete_model(&d4_model);
    d4_model = disprs_d4_new_model(d4_error, d4_structure, DISPRS_D4S);
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);
    assert(!disprs_d4_check_error(d4_error));
    assert(isfinite(energy) && energy < 0.0);

    disprs_d4_delete_param(&d4_param);
    disprs_d4_delete_model(&d4_model);
    disprs_d4_delete_structure(&d4_structure);
    disprs_d4_delete_error(&d4_error);

    d4_error = disprs_d4_new_error();
    d4_structure = disprs_d4_new_structure(d4_error, 2, periodic_numbers, periodic_positions, NULL,
                                           lattice, periodic);
    d4_model = disprs_d4_new_model(d4_error, d4_structure, DISPRS_D4S);
    d4_param = disprs_d4_load_param(d4_error, "pbe", false);
    disprs_d4_set_model_ewald(d4_error, d4_model, 0, 1e-10, 0.0, -1);
    assert(!disprs_d4_check_error(d4_error));
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);
    assert(!disprs_d4_check_error(d4_error) && isfinite(energy));
    const double direct_energy = energy;
    disprs_d4_set_model_ewald(d4_error, d4_model, 0, NAN, 0.0, 32);
    assert(disprs_d4_check_error(d4_error));
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, NULL, NULL);
    assert(!disprs_d4_check_error(d4_error) && energy == direct_energy);
    disprs_d4_set_model_ewald(d4_error, d4_model, 0, 1e-10, 0.0, 32);
    disprs_d4_get_dispersion(d4_error, d4_structure, d4_model, d4_param, &energy, gradient, virial);
    assert(!disprs_d4_check_error(d4_error) && fabs(energy - direct_energy) < 5e-8);
    analytical_hessian[0] = numerical_hessian[0] = pair2[0] = 42.0;
    disprs_d4_get_dispersion_hessian(d4_error, d4_structure, d4_model, d4_param,
                                     analytical_hessian);
    assert(disprs_d4_check_error(d4_error) && analytical_hessian[0] == 42.0);
    disprs_d4_get_numerical_hessian(d4_error, d4_structure, d4_model, d4_param, numerical_hessian);
    assert(disprs_d4_check_error(d4_error) && numerical_hessian[0] == 42.0);
    disprs_d4_get_pairwise_dispersion(d4_error, d4_structure, d4_model, d4_param, pair2, pair3);
    assert(disprs_d4_check_error(d4_error) && pair2[0] == 42.0);
    disprs_d4_delete_param(&d4_param);
    disprs_d4_delete_model(&d4_model);
    disprs_d4_delete_structure(&d4_structure);
    disprs_d4_delete_error(&d4_error);
    return 0;
}
