/**************************************************************************

   Copyright © 2018 Sergey Poluyan <svpoluyan@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

**************************************************************************/
#ifndef INCLUDED_Dunbrackdata_hh
#define INCLUDED_Dunbrackdata_hh

#include <string>
#include <vector>
#include <iostream>
#include <random>

namespace pepsgo
{
	namespace bbdep
	{

/// designation
// 60 180 300
// 60 180 -60
// 1  2   3
// gauche+ (g+), trans (t), and gauche- (g-)

/// Rosetta
//0 AG
//1	V
//2 WHNDFILCTS
//3 MEQPY
//4 RK

/// Dunbrack [Number of bins for each discrete chi angle]
//1 SVCT [3] ser val cys thr
//2 LIFDYWNH [3, 3] leu ile    [3, 6] phe asp tyr     [3, 12] trp asn his
//3 MEQP [3, 3, 3] met          [3, 3, 6] glu         [3, 3, 12] gln           [2, 1, 1] pro
//4 RK [3, 3, 3, 3, 1] arg [3, 3, 3, 3] lys

/// Dunbrack rotameric
//	ARG χ1, χ2, χ3, χ4    HIS χ1                PRO χ1
//	ASN χ1                ILE χ1                SER χ1
//	ASP χ1                LEU χ1, χ2            THR χ1
//	CYS χ1                LYS χ1, χ2, χ3, χ4    TRP χ1
//	GLN χ1, χ2            MET χ1, χ2, χ3        TYR χ1
//	GLU χ1, χ2            PHE χ1                VAL χ1

/// Dunbrack non-rotameric
//	ASN χ2    GLN χ3    PHE χ2    HIS χ2
//	ASP χ2    GLU χ3    TYR χ2    TRP χ2

		struct Dunbrack_data
		{
			std::string name;
			double Phi;
			double Psi;
			int Count;
			size_t r1;
			size_t r2;
			size_t r3;
			size_t r4;
			double Probabil;
			double chi1Val;
			double chi2Val;
			double chi3Val;
			double chi4Val;

			double chi1Sig;
			double chi2Sig;
			double chi3Sig;
			double chi4Sig;

			friend std::ostream& operator<<(std::ostream& stream, const Dunbrack_data& obj);
		};

		std::string remove_chars(const std::string is, const std::string& chars);

		void load_data_sm(std::string fname, std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers);

		Dunbrack_data get_max_prob_object(std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers, double Phi, double Psi);

		Dunbrack_data get_probabil(std::vector<Dunbrack_data> &data, std::vector<std::vector<int>> &numbers, std::mt19937_64 &generator, double Phi, double Psi);

		bool has_only_zeros_or_dots(const std::string str);

	}
}

#endif
