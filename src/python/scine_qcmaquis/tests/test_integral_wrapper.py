import os

import numpy as np
import pytest

from scine_qcmaquis.integral_wrapper import IntegralMapWrapper, IntegralsParser, IntegralUtils

class TestIntegralWrapper:
    def test_integral_utils_permute_same_particle(self):
        int_utils = IntegralUtils()
        result_list = [[0, 1, 2, 3], [0, 1, 3, 2]]
        result_list = int_utils._permute_same_particle(result_list, 0, 1)
        test_list = [[0, 1, 2, 3], [0, 1, 3, 2], [1, 0, 2, 3], [1, 0, 3, 2]]
        assert result_list == test_list


    def test_integral_utils_permute_particle_block(self):
        int_utils = IntegralUtils()
        result_list = [[0, 1, 2, 3], [0, 1, 3, 2]]
        result_list = int_utils._permute_particle_block(result_list, 0, 2)
        test_list = [[0, 1, 2, 3], [0, 1, 3, 2], [2, 3, 0, 1], [3, 2, 0, 1]]
        assert result_list == test_list


    def test_integral_utils_get_symmetric_indices(self):
        int_utils = IntegralUtils()
        result_list_1 = [0, 1, 2, 3]
        result_list_1 = int_utils.get_symmetric_indices(result_list_1, "one")
        test_list_1 = [[0, 1, 2, 3]]
        assert result_list_1 == test_list_1

        result_list_2 = [0, 1, 2, 3]
        result_list_2 = int_utils.get_symmetric_indices(result_list_2, "two")
        test_list_2 = [[0, 1, 2, 3], [2, 3, 0, 1]]
        assert result_list_2 == test_list_2

        result_list_4 = [0, 1, 2, 3]
        result_list_4 = int_utils.get_symmetric_indices(result_list_4, "four")
        test_list_4 = [[0, 1, 2, 3], [1, 0, 2, 3], [0, 1, 3, 2], [1, 0, 3, 2]]
        assert result_list_4 == test_list_4

        result_list_8 = [0, 1, 2, 3]
        result_list_8 = int_utils.get_symmetric_indices(result_list_8, "eight")
        test_list_8 = [[0, 1, 2, 3], [2, 3, 0, 1], [1, 0, 2, 3], [3, 2, 0, 1], [0, 1, 3, 2], [2, 3, 1, 0], [1, 0, 3, 2], [3, 2, 1, 0]]
        assert result_list_8 == test_list_8


    def test_integral_parser_parse_fcidump_header(self):
        def _write_fcidump():
            tmp_fcidump = open("fcidump_mock", "w")
            tmp_fcidump.write("&FCI NORB=19, NELEC=4, MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.close()
        _write_fcidump()

        parser = IntegralsParser()
        fcidump = open("fcidump_mock", "r")
        parser._parse_fcidump_header(fcidump)
        fcidump.close()
        os.remove("fcidump_mock")

        assert parser.fcidump_values.norb == 19
        assert parser.fcidump_values.nelec == 4
        assert parser.fcidump_values.ms2 == 0
        assert parser.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert parser.fcidump_values.isym == 1
        assert parser.fcidump_values.transcorrelated is False
        assert parser.fcidump_values.unrestricted is False

        def _write_fcidump_2():
            tmp_fcidump = open("fcidump_mock_2", "w")
            tmp_fcidump.write("&FCI\n")
            tmp_fcidump.write("NORB=19,\n")
            tmp_fcidump.write("NELEC=4,\n")
            tmp_fcidump.write("MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.close()
        _write_fcidump_2()

        parser_2 = IntegralsParser()
        fcidump_2 = open("fcidump_mock_2", "r")
        parser_2._parse_fcidump_header(fcidump_2)
        fcidump_2.close()
        os.remove("fcidump_mock_2")

        assert parser_2.fcidump_values.norb == 19
        assert parser_2.fcidump_values.nelec == 4
        assert parser_2.fcidump_values.ms2 == 0
        assert parser_2.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert parser_2.fcidump_values.isym == 1
        assert parser_2.fcidump_values.transcorrelated is False
        assert parser_2.fcidump_values.unrestricted is False


    def test_integral_parser_parse_fcidump_body_incorrect_symmetry_in_fcidump(self):
        def _write_fcidump():
            tmp_fcidump = open("fcidump_mock", "w")
            tmp_fcidump.write("&FCI NORB=19, NELEC=4, MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.write("  -2.4500540021458295e+00   1   1   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   1   2   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   2   1   0   0\n")
            tmp_fcidump.write("   3.6933459762744647e-12   1   3   0   0\n")
            tmp_fcidump.write("  -1.1965308951977782e-11   1   6   0   0\n")
            tmp_fcidump.write("   1.6482623985379075e+00   1   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   2   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   1   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   2   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   1   2   1   2\n")
            tmp_fcidump.write("   3.4283504997827774e-01   2   2   1   1\n")
            tmp_fcidump.write("   6.4104605406089948e-03   2   2   2   1\n")
            tmp_fcidump.write("   4.7954845474515156e-01   2   2   2   2\n")
            tmp_fcidump.write("  -6.0591020065061930e-02   3   1   1   1\n")
            tmp_fcidump.write("   4.2163595616380176e-03   3   1   2   1\n")
            tmp_fcidump.write("  -5.7500401969882015e-03   3   1   2   2\n")
            tmp_fcidump.write("   4.2203448717115129e-03   3   1   3   1\n")
            tmp_fcidump.write("   3.6452515524633749e-03   3   2   1   1\n")
            tmp_fcidump.write("  -1.7606280599412798e-03   3   2   2   1\n")
            tmp_fcidump.write("  -4.0131224007160043e-02   3   2   2   2\n")
            tmp_fcidump.write("   4.0416762395947006e-04   3   2   3   1\n")
            tmp_fcidump.write("   7.3361564717230529e-03   3   2   3   2\n")
            tmp_fcidump.write("   2.3122843673950771e-01   3   3   1   1\n")
            tmp_fcidump.write("  -3.0256008119892365e-03   3   3   2   1\n")
            tmp_fcidump.write("   1.5763733543509206e-01   3   3   2   2\n")
            tmp_fcidump.write("   1.3138548086074515e-03   3   3   3   1\n")
            tmp_fcidump.write("   6.1029726451229879e-03   3   3   3   2\n")
            tmp_fcidump.write("   2.0527241559476508e-01   3   3   3   3\n")
            tmp_fcidump.write("   2.7912933884382570e-04   4   1   4   1\n")
            tmp_fcidump.write("   4.8281756725673466e-04   4   2   4   1\n")
            tmp_fcidump.write("   5.4383779920285152e-03   4   2   4   2\n")
            tmp_fcidump.write("   7.5847612146573342e-04   4   3   4   1\n")
            tmp_fcidump.write("   6.4842523433864055e-03   4   3   4   2\n")
            tmp_fcidump.write("   2.1099942529175689e-02   4   3   4   3\n")
            tmp_fcidump.write("   1.7497909517853510e-01   4   4   1   1\n")
            tmp_fcidump.write("  -3.7240137619155936e-04   4   4   2   1\n")
            tmp_fcidump.write("   1.5589038649219863e-01   4   4   2   2\n")
            tmp_fcidump.write("   1.2949829141455018e-04   4   4   3   1\n")
            tmp_fcidump.write("   1.8750572872915801e-03   4   4   3   2\n")
            tmp_fcidump.write("   1.4661471245002852e-01   4   4   3   3\n")
            tmp_fcidump.write("   1.4153931193989575e-01   4   4   4   4\n")
            tmp_fcidump.close()
        _write_fcidump()

        parser = IntegralsParser()
        parser.parse_fcidump("fcidump_mock")

        assert parser.fcidump_values.norb == 19
        assert parser.fcidump_values.nelec == 4
        assert parser.fcidump_values.ms2 == 0
        assert parser.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert parser.fcidump_values.isym == 1
        assert parser.fcidump_values.transcorrelated is False
        assert parser.fcidump_values.unrestricted is False
        unique_terms = {
            (0, 0, 0, 0): -7.9837309384386277e+00,
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-12,
            (1, 3, 0, 0): 3.6933459762744647e-12,
            (1, 6, 0, 0): -1.1965308951977782e-11,
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (2, 1, 1, 1): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        assert unique_terms == parser._unique_term


    def test_integral_parser_parse_fcidump_body(self):
        def _write_fcidump():
            tmp_fcidump = open("fcidump_mock", "w")
            tmp_fcidump.write("&FCI NORB=19, NELEC=4, MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.write("  -2.4500540021458295e+00   1   1   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-02   1   2   0   0\n")
            tmp_fcidump.write("   3.6933459762744647e-02   1   3   0   0\n")
            tmp_fcidump.write("  -1.1965308951977782e-01   1   6   0   0\n")
            tmp_fcidump.write("   1.6482623985379075e+00   1   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   2   1   1   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   2   1   2   1\n")
            tmp_fcidump.write("   3.4283504997827774e-01   2   2   1   1\n")
            tmp_fcidump.write("   6.4104605406089948e-03   2   2   2   1\n")
            tmp_fcidump.write("   4.7954845474515156e-01   2   2   2   2\n")
            tmp_fcidump.write("  -6.0591020065061930e-02   3   1   1   1\n")
            tmp_fcidump.write("   4.2163595616380176e-03   3   1   2   1\n")
            tmp_fcidump.write("  -5.7500401969882015e-03   3   1   2   2\n")
            tmp_fcidump.write("   4.2203448717115129e-03   3   1   3   1\n")
            tmp_fcidump.write("   3.6452515524633749e-03   3   2   1   1\n")
            tmp_fcidump.write("  -1.7606280599412798e-03   3   2   2   1\n")
            tmp_fcidump.write("  -4.0131224007160043e-02   3   2   2   2\n")
            tmp_fcidump.write("   4.0416762395947006e-04   3   2   3   1\n")
            tmp_fcidump.write("   7.3361564717230529e-03   3   2   3   2\n")
            tmp_fcidump.write("   2.3122843673950771e-01   3   3   1   1\n")
            tmp_fcidump.write("  -3.0256008119892365e-03   3   3   2   1\n")
            tmp_fcidump.write("   1.5763733543509206e-01   3   3   2   2\n")
            tmp_fcidump.write("   1.3138548086074515e-03   3   3   3   1\n")
            tmp_fcidump.write("   6.1029726451229879e-03   3   3   3   2\n")
            tmp_fcidump.write("   2.0527241559476508e-01   3   3   3   3\n")
            tmp_fcidump.write("   2.7912933884382570e-04   4   1   4   1\n")
            tmp_fcidump.write("   4.8281756725673466e-04   4   2   4   1\n")
            tmp_fcidump.write("   5.4383779920285152e-03   4   2   4   2\n")
            tmp_fcidump.write("   7.5847612146573342e-04   4   3   4   1\n")
            tmp_fcidump.write("   6.4842523433864055e-03   4   3   4   2\n")
            tmp_fcidump.write("   2.1099942529175689e-02   4   3   4   3\n")
            tmp_fcidump.write("   1.7497909517853510e-01   4   4   1   1\n")
            tmp_fcidump.write("  -3.7240137619155936e-04   4   4   2   1\n")
            tmp_fcidump.write("   1.5589038649219863e-01   4   4   2   2\n")
            tmp_fcidump.write("   1.2949829141455018e-04   4   4   3   1\n")
            tmp_fcidump.write("   1.8750572872915801e-03   4   4   3   2\n")
            tmp_fcidump.write("   1.4661471245002852e-01   4   4   3   3\n")
            tmp_fcidump.write("   1.4153931193989575e-01   4   4   4   4\n")
            tmp_fcidump.close()
        _write_fcidump()

        parser = IntegralsParser()
        parser.parse_fcidump("fcidump_mock")

        assert parser.fcidump_values.norb == 19
        assert parser.fcidump_values.nelec == 4
        assert parser.fcidump_values.ms2 == 0
        assert parser.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert parser.fcidump_values.isym == 1
        assert parser.fcidump_values.transcorrelated is False
        assert parser.fcidump_values.unrestricted is False
        unique_terms = {
            (0, 0, 0, 0): -7.9837309384386277e+00,
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-02,
            (1, 3, 0, 0): 3.6933459762744647e-02,
            (1, 6, 0, 0): -1.1965308951977782e-01,
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (2, 1, 1, 1): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        assert unique_terms == parser._unique_term


    def test_integral_set_core(self):
        core_value = -0.983267904
        unique_terms = {(0, 0, 0, 0): core_value}
        parser = IntegralsParser()
        parser.set_core(core_value)
        assert unique_terms == parser._unique_term


    def test_integral_parse_one_body(self):
        one_body = np.zeros((6, 6))
        one_body[0, 0] = -2.4500540021458295e+00
        one_body[0, 1] = -7.2269742273922688e-02
        one_body[1, 0] = -7.2269742273922688e-02
        one_body[0, 2] = 3.6933459762744647e-02
        one_body[2, 0] = 3.6933459762744647e-02
        one_body[0, 5] = -1.1965308951977782e-01
        one_body[5, 0] = -1.1965308951977782e-01

        unique_terms = {
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-02,
            (1, 3, 0, 0): 3.6933459762744647e-02,
            (1, 6, 0, 0): -1.1965308951977782e-01,
        }
        parser = IntegralsParser()
        parser.parse_one_body(one_body, 6)
        assert unique_terms == parser._unique_term


    def test_integral_parse_two_body(self):
        two_body = np.zeros((4, 4, 4, 4))
        two_body[0, 0, 0, 0] = 1.6482623985379075e+00
        two_body[1, 0, 0, 0] = -8.9231283142061496e-02
        two_body[0, 1, 0, 0] = -8.9231283142061496e-02
        two_body[0, 0, 1, 0] = -8.9231283142061496e-02
        two_body[0, 0, 0, 1] = -8.9231283142061496e-02
        two_body[1, 0, 1, 0] = 8.7306915386681453e-03
        two_body[1, 1, 0, 0] = 3.4283504997827774e-01
        two_body[1, 1, 1, 0] = 6.4104605406089948e-03
        two_body[1, 1, 1, 1] = 4.7954845474515156e-01
        two_body[2, 0, 0, 0] = -6.0591020065061930e-02
        two_body[2, 0, 1, 0] = 4.2163595616380176e-03
        two_body[2, 0, 1, 1] = -5.7500401969882015e-03
        two_body[2, 0, 2, 0] = 4.2203448717115129e-03
        two_body[2, 1, 0, 0] = 3.6452515524633749e-03
        two_body[2, 1, 1, 0] = -1.7606280599412798e-03
        two_body[2, 1, 1, 1] = -4.0131224007160043e-02
        two_body[2, 1, 2, 0] = 4.0416762395947006e-04
        two_body[2, 1, 2, 1] = 7.3361564717230529e-03
        two_body[2, 2, 0, 0] = 2.3122843673950771e-01
        two_body[2, 2, 1, 0] = -3.0256008119892365e-03
        two_body[2, 2, 1, 1] = 1.5763733543509206e-01
        two_body[2, 2, 2, 0] = 1.3138548086074515e-03
        two_body[2, 2, 2, 1] = 6.1029726451229879e-03
        two_body[2, 2, 2, 2] = 2.0527241559476508e-01
        two_body[3, 0, 3, 0] = 2.7912933884382570e-04
        two_body[3, 1, 3, 0] = 4.8281756725673466e-04
        two_body[3, 1, 3, 1] = 5.4383779920285152e-03
        two_body[3, 2, 3, 0] = 7.5847612146573342e-04
        two_body[3, 2, 3, 1] = 6.4842523433864055e-03
        two_body[3, 2, 3, 2] = 2.1099942529175689e-02
        two_body[3, 3, 0, 0] = 1.7497909517853510e-01
        two_body[3, 3, 1, 0] = -3.7240137619155936e-04
        two_body[3, 3, 1, 1] = 1.5589038649219863e-01
        two_body[3, 3, 2, 0] = 1.2949829141455018e-04
        two_body[3, 3, 2, 1] = 1.8750572872915801e-03
        two_body[3, 3, 2, 2] = 1.4661471245002852e-01
        two_body[3, 3, 3, 3] = 1.4153931193989575e-01

        unique_terms = {
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (1, 1, 1, 2): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        parser = IntegralsParser()
        parser.parse_two_body(two_body, 4)
        assert unique_terms == parser._unique_term


    def test_integral_map_wrapper_pyscf(self):
        core_value = -0.983267904

        one_body = np.zeros((4, 4))
        one_body[0, 0] = -2.4500540021458295e+00
        one_body[0, 1] = -7.2269742273922688e-02
        one_body[1, 0] = -7.2269742273922688e-02
        one_body[0, 2] = 3.6933459762744647e-02
        one_body[2, 0] = 3.6933459762744647e-02
        one_body[2, 3] = -1.1965308951977782e-01
        one_body[3, 2] = -1.1965308951977782e-01

        two_body = np.zeros((4, 4, 4, 4))
        two_body[0, 0, 0, 0] = 1.6482623985379075e+00
        two_body[1, 0, 0, 0] = -8.9231283142061496e-02
        two_body[0, 1, 0, 0] = -8.9231283142061496e-02
        two_body[0, 0, 1, 0] = -8.9231283142061496e-02
        two_body[0, 0, 0, 1] = -8.9231283142061496e-02
        two_body[1, 0, 1, 0] = 8.7306915386681453e-03
        two_body[1, 1, 0, 0] = 3.4283504997827774e-01
        two_body[1, 1, 1, 0] = 6.4104605406089948e-03
        two_body[1, 1, 1, 1] = 4.7954845474515156e-01
        two_body[2, 0, 0, 0] = -6.0591020065061930e-02
        two_body[2, 0, 1, 0] = 4.2163595616380176e-03
        two_body[2, 0, 1, 1] = -5.7500401969882015e-03
        two_body[2, 0, 2, 0] = 4.2203448717115129e-03
        two_body[2, 1, 0, 0] = 3.6452515524633749e-03
        two_body[2, 1, 1, 0] = -1.7606280599412798e-03
        two_body[2, 1, 1, 1] = -4.0131224007160043e-02
        two_body[2, 1, 2, 0] = 4.0416762395947006e-04
        two_body[2, 1, 2, 1] = 7.3361564717230529e-03
        two_body[2, 2, 0, 0] = 2.3122843673950771e-01
        two_body[2, 2, 1, 0] = -3.0256008119892365e-03
        two_body[2, 2, 1, 1] = 1.5763733543509206e-01
        two_body[2, 2, 2, 0] = 1.3138548086074515e-03
        two_body[2, 2, 2, 1] = 6.1029726451229879e-03
        two_body[2, 2, 2, 2] = 2.0527241559476508e-01
        two_body[3, 0, 3, 0] = 2.7912933884382570e-04
        two_body[3, 1, 3, 0] = 4.8281756725673466e-04
        two_body[3, 1, 3, 1] = 5.4383779920285152e-03
        two_body[3, 2, 3, 0] = 7.5847612146573342e-04
        two_body[3, 2, 3, 1] = 6.4842523433864055e-03
        two_body[3, 2, 3, 2] = 2.1099942529175689e-02
        two_body[3, 3, 0, 0] = 1.7497909517853510e-01
        two_body[3, 3, 1, 0] = -3.7240137619155936e-04
        two_body[3, 3, 1, 1] = 1.5589038649219863e-01
        two_body[3, 3, 2, 0] = 1.2949829141455018e-04
        two_body[3, 3, 2, 1] = 1.8750572872915801e-03
        two_body[3, 3, 2, 2] = 1.4661471245002852e-01
        two_body[3, 3, 3, 3] = 1.4153931193989575e-01

        unique_terms = {
            (0, 0, 0, 0): core_value,
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-02,
            (1, 3, 0, 0): 3.6933459762744647e-02,
            (3, 4, 0, 0): -1.1965308951977782e-01,
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (1, 1, 1, 2): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        map_wrapper = IntegralMapWrapper()
        norb = 4
        map_wrapper.fill_from_pyscf(core_value, one_body, two_body, norb)
        assert unique_terms == map_wrapper._parser._unique_term


    def test_integral_wrapper_fcidump_body_incorrect_symmetry_in_fcidump(self):
        def _write_fcidump():
            tmp_fcidump = open("fcidump_mock", "w")
            tmp_fcidump.write("&FCI NORB=19, NELEC=4, MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.write("  -2.4500540021458295e+00   1   1   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   1   2   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   2   1   0   0\n")
            tmp_fcidump.write("   3.6933459762744647e-12   1   3   0   0\n")
            tmp_fcidump.write("  -1.1965308951977782e-11   1   6   0   0\n")
            tmp_fcidump.write("   1.6482623985379075e+00   1   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   2   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   1   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   2   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   1   2   1   2\n")
            tmp_fcidump.write("   3.4283504997827774e-01   2   2   1   1\n")
            tmp_fcidump.write("   6.4104605406089948e-03   2   2   2   1\n")
            tmp_fcidump.write("   4.7954845474515156e-01   2   2   2   2\n")
            tmp_fcidump.write("  -6.0591020065061930e-02   3   1   1   1\n")
            tmp_fcidump.write("   4.2163595616380176e-03   3   1   2   1\n")
            tmp_fcidump.write("  -5.7500401969882015e-03   3   1   2   2\n")
            tmp_fcidump.write("   4.2203448717115129e-03   3   1   3   1\n")
            tmp_fcidump.write("   3.6452515524633749e-03   3   2   1   1\n")
            tmp_fcidump.write("  -1.7606280599412798e-03   3   2   2   1\n")
            tmp_fcidump.write("  -4.0131224007160043e-02   3   2   2   2\n")
            tmp_fcidump.write("   4.0416762395947006e-04   3   2   3   1\n")
            tmp_fcidump.write("   7.3361564717230529e-03   3   2   3   2\n")
            tmp_fcidump.write("   2.3122843673950771e-01   3   3   1   1\n")
            tmp_fcidump.write("  -3.0256008119892365e-03   3   3   2   1\n")
            tmp_fcidump.write("   1.5763733543509206e-01   3   3   2   2\n")
            tmp_fcidump.write("   1.3138548086074515e-03   3   3   3   1\n")
            tmp_fcidump.write("   6.1029726451229879e-03   3   3   3   2\n")
            tmp_fcidump.write("   2.0527241559476508e-01   3   3   3   3\n")
            tmp_fcidump.write("   2.7912933884382570e-04   4   1   4   1\n")
            tmp_fcidump.write("   4.8281756725673466e-04   4   2   4   1\n")
            tmp_fcidump.write("   5.4383779920285152e-03   4   2   4   2\n")
            tmp_fcidump.write("   7.5847612146573342e-04   4   3   4   1\n")
            tmp_fcidump.write("   6.4842523433864055e-03   4   3   4   2\n")
            tmp_fcidump.write("   2.1099942529175689e-02   4   3   4   3\n")
            tmp_fcidump.write("   1.7497909517853510e-01   4   4   1   1\n")
            tmp_fcidump.write("  -3.7240137619155936e-04   4   4   2   1\n")
            tmp_fcidump.write("   1.5589038649219863e-01   4   4   2   2\n")
            tmp_fcidump.write("   1.2949829141455018e-04   4   4   3   1\n")
            tmp_fcidump.write("   1.8750572872915801e-03   4   4   3   2\n")
            tmp_fcidump.write("   1.4661471245002852e-01   4   4   3   3\n")
            tmp_fcidump.write("   1.4153931193989575e-01   4   4   4   4\n")
            tmp_fcidump.close()
        _write_fcidump()

        parser = IntegralsParser()
        parser.parse_fcidump("fcidump_mock")

        assert parser.fcidump_values.norb == 19
        assert parser.fcidump_values.nelec == 4
        assert parser.fcidump_values.ms2 == 0
        assert parser.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert parser.fcidump_values.isym == 1
        assert parser.fcidump_values.transcorrelated is False
        assert parser.fcidump_values.unrestricted is False
        unique_terms = {
            (0, 0, 0, 0): -7.9837309384386277e+00,
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-12,
            (1, 3, 0, 0): 3.6933459762744647e-12,
            (1, 6, 0, 0): -1.1965308951977782e-11,
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (2, 1, 1, 1): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        assert unique_terms == parser._unique_term


    def test_integral_wrapper_fcidump(self):
        def _write_fcidump():
            tmp_fcidump = open("fcidump_mock", "w")
            tmp_fcidump.write("&FCI NORB=19, NELEC=4, MS2=0,\n")
            tmp_fcidump.write("ORBSYM=1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,\n")
            tmp_fcidump.write("ISYM=1,\n")
            tmp_fcidump.write("&END\n")
            tmp_fcidump.write("  -7.9837309384386277e+00   0   0   0   0\n")
            tmp_fcidump.write("  -2.4500540021458295e+00   1   1   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   1   2   0   0\n")
            tmp_fcidump.write("  -7.2269742273922688e-12   2   1   0   0\n")
            tmp_fcidump.write("   3.6933459762744647e-12   1   3   0   0\n")
            tmp_fcidump.write("  -1.1965308951977782e-11   1   6   0   0\n")
            tmp_fcidump.write("   1.6482623985379075e+00   1   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   2   1   1   1\n")
            tmp_fcidump.write("  -8.9231283142061496e-02   1   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   2   1   2   1\n")
            tmp_fcidump.write("   8.7306915386681453e-03   1   2   1   2\n")
            tmp_fcidump.write("   3.4283504997827774e-01   2   2   1   1\n")
            tmp_fcidump.write("   6.4104605406089948e-03   2   2   2   1\n")
            tmp_fcidump.write("   4.7954845474515156e-01   2   2   2   2\n")
            tmp_fcidump.write("  -6.0591020065061930e-02   3   1   1   1\n")
            tmp_fcidump.write("   4.2163595616380176e-03   3   1   2   1\n")
            tmp_fcidump.write("  -5.7500401969882015e-03   3   1   2   2\n")
            tmp_fcidump.write("   4.2203448717115129e-03   3   1   3   1\n")
            tmp_fcidump.write("   3.6452515524633749e-03   3   2   1   1\n")
            tmp_fcidump.write("  -1.7606280599412798e-03   3   2   2   1\n")
            tmp_fcidump.write("  -4.0131224007160043e-02   3   2   2   2\n")
            tmp_fcidump.write("   4.0416762395947006e-04   3   2   3   1\n")
            tmp_fcidump.write("   7.3361564717230529e-03   3   2   3   2\n")
            tmp_fcidump.write("   2.3122843673950771e-01   3   3   1   1\n")
            tmp_fcidump.write("  -3.0256008119892365e-03   3   3   2   1\n")
            tmp_fcidump.write("   1.5763733543509206e-01   3   3   2   2\n")
            tmp_fcidump.write("   1.3138548086074515e-03   3   3   3   1\n")
            tmp_fcidump.write("   6.1029726451229879e-03   3   3   3   2\n")
            tmp_fcidump.write("   2.0527241559476508e-01   3   3   3   3\n")
            tmp_fcidump.write("   2.7912933884382570e-04   4   1   4   1\n")
            tmp_fcidump.write("   4.8281756725673466e-04   4   2   4   1\n")
            tmp_fcidump.write("   5.4383779920285152e-03   4   2   4   2\n")
            tmp_fcidump.write("   7.5847612146573342e-04   4   3   4   1\n")
            tmp_fcidump.write("   6.4842523433864055e-03   4   3   4   2\n")
            tmp_fcidump.write("   2.1099942529175689e-02   4   3   4   3\n")
            tmp_fcidump.write("   1.7497909517853510e-01   4   4   1   1\n")
            tmp_fcidump.write("  -3.7240137619155936e-04   4   4   2   1\n")
            tmp_fcidump.write("   1.5589038649219863e-01   4   4   2   2\n")
            tmp_fcidump.write("   1.2949829141455018e-04   4   4   3   1\n")
            tmp_fcidump.write("   1.8750572872915801e-03   4   4   3   2\n")
            tmp_fcidump.write("   1.4661471245002852e-01   4   4   3   3\n")
            tmp_fcidump.write("   1.4153931193989575e-01   4   4   4   4\n")
            tmp_fcidump.close()
        _write_fcidump()

        integral_map = IntegralMapWrapper()
        integral_map.fill_from_fcidump("fcidump_mock")

        assert integral_map._parser.fcidump_values.norb == 19
        assert integral_map._parser.fcidump_values.nelec == 4
        assert integral_map._parser.fcidump_values.ms2 == 0
        assert integral_map._parser.fcidump_values.orbsym == [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        assert integral_map._parser.fcidump_values.isym == 1
        assert integral_map._parser.fcidump_values.transcorrelated is False
        assert integral_map._parser.fcidump_values.unrestricted is False
        unique_terms = {
            (0, 0, 0, 0): -7.9837309384386277e+00,
            (1, 1, 0, 0): -2.4500540021458295e+00,
            (1, 2, 0, 0): -7.2269742273922688e-12,
            (1, 3, 0, 0): 3.6933459762744647e-12,
            (1, 6, 0, 0): -1.1965308951977782e-11,
            (1, 1, 1, 1): 1.6482623985379075e+00,
            (2, 1, 1, 1): -8.9231283142061496e-02,
            (2, 1, 2, 1): 8.7306915386681453e-03,
            (2, 2, 1, 1): 3.4283504997827774e-01,
            (2, 2, 2, 1): 6.4104605406089948e-03,
            (2, 2, 2, 2): 4.7954845474515156e-01,
            (3, 1, 1, 1): -6.0591020065061930e-02,
            (3, 1, 2, 1): 4.2163595616380176e-03,
            (3, 1, 2, 2): -5.7500401969882015e-03,
            (3, 1, 3, 1): 4.2203448717115129e-03,
            (3, 2, 1, 1): 3.6452515524633749e-03,
            (3, 2, 2, 1): -1.7606280599412798e-03,
            (3, 2, 2, 2): -4.0131224007160043e-02,
            (3, 2, 3, 1): 4.0416762395947006e-04,
            (3, 2, 3, 2): 7.3361564717230529e-03,
            (3, 3, 1, 1): 2.3122843673950771e-01,
            (3, 3, 2, 1): -3.0256008119892365e-03,
            (3, 3, 2, 2): 1.5763733543509206e-01,
            (3, 3, 3, 1): 1.3138548086074515e-03,
            (3, 3, 3, 2): 6.1029726451229879e-03,
            (3, 3, 3, 3): 2.0527241559476508e-01,
            (4, 1, 4, 1): 2.7912933884382570e-04,
            (4, 2, 4, 1): 4.8281756725673466e-04,
            (4, 2, 4, 2): 5.4383779920285152e-03,
            (4, 3, 4, 1): 7.5847612146573342e-04,
            (4, 3, 4, 2): 6.4842523433864055e-03,
            (4, 3, 4, 3): 2.1099942529175689e-02,
            (4, 4, 1, 1): 1.7497909517853510e-01,
            (4, 4, 2, 1): -3.7240137619155936e-04,
            (4, 4, 2, 2): 1.5589038649219863e-01,
            (4, 4, 3, 1): 1.2949829141455018e-04,
            (4, 4, 3, 2): 1.8750572872915801e-03,
            (4, 4, 3, 3): 1.4661471245002852e-01,
            (4, 4, 4, 4): 1.4153931193989575e-01,
        }
        assert unique_terms == integral_map._parser._unique_term


if __name__ == "__main__":
    test_integral_utils_permute_same_particle()
    test_integral_utils_permute_particle_block()
    test_integral_utils_get_symmetric_indices()
    test_integral_parser_parse_fcidump_header()
    test_integral_parser_parse_fcidump_body()
    test_integral_parser_parse_fcidump_body_incorrect_symmetry_in_fcidump()
    test_integral_set_core()
    test_integral_parse_one_body()
    test_integral_parse_two_body()
    test_integral_map_wrapper_pyscf()
    test_integral_wrapper_fcidump()
