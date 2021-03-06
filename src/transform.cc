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
#include <transform.hh>
#include <bbtools.hh>

#include <numeric/conversions.hh>

namespace pepsgo
{
	namespace transform
	{
		void bbdep_experiment_actual_states(
		  const std::vector<double> &x,
		  std::vector<double> &out,
		  const std::vector<pepsgo::opt_element> &opt_vect,
		  const pepsgo::ranges &range,
		  const pepsgo::bbdep::BBDEP_Dunbrack_sm &bbdep_obj_sm,
		  size_t peptide_phipsi_2d_size)
		{
			if(!range.do_chi) // bbdep_experiment
				return;

			if(std::get<0>(range.chi) == 0)
				return;

			for(size_t i = std::get<1>(range.chi); i <= std::get<2>(range.chi); i++)
			{
				out[i] = -numeric::NumericTraits<core::Real>::pi()*(1.0 - 2.0*x[i]);
			}
			double phi = 0, psi = 0;
			for(size_t c = std::get<1>(range.chi); c <= std::get<2>(range.chi);)
			{
				pepsgo::bbdep::Dunbrack_data first_line = bbdep_obj_sm.get_first_line(opt_vect[c].amino_acid);

				if(opt_vect[c].seqpos == 1)
				{
					if(opt_vect[c].nchi == 1)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
					}
					else if(opt_vect[c].nchi == 2)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}
					}
					else if(opt_vect[c].nchi == 3)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}

						if(first_line.r3 != 0)
						{
							double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 2], 2);
							out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
						}
					}
					else if(opt_vect[c].nchi == 4)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}

						if(first_line.r3 != 0)
						{
							double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 2], 2);
							out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
						}

						if(first_line.r4 != 0)
						{
							double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 3], 3);
							out[c + 3] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
						}
					}

					c += opt_vect[c].nchi;
					continue;
				}
				else if(opt_vect[c].seqpos == peptide_phipsi_2d_size + 2)
				{
					if(opt_vect[c].nchi == 1)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
					}
					else if(opt_vect[c].nchi == 2)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}
					}
					else if(opt_vect[c].nchi == 3)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}

						if(first_line.r3 != 0)
						{
							double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 2], 2);
							out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
						}
					}
					else if(opt_vect[c].nchi == 4)
					{
						double c1 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c], 0);
						out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

						if(first_line.r2 != 0)
						{
							double c2 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 1], 1);
							out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
						}

						if(first_line.r3 != 0)
						{
							double c3 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 2], 2);
							out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
						}

						if(first_line.r4 != 0)
						{
							double c4 = bbdep_obj_sm.get_degree_bbind(opt_vect[c].amino_acid, x[c + 3], 3);
							out[c + 3] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
						}
					}

					c += opt_vect[c].nchi;
					continue;
				}
				else
				{
					core::Size sp = opt_vect[c].seqpos;
					auto it = std::find_if(opt_vect.begin(), opt_vect.begin() + 2 * peptide_phipsi_2d_size + 2, [&sp](const opt_element &s)
					{
						return s.seqpos == sp && s.torsion_name == "bbphi";
					});
					size_t index = std::distance(opt_vect.begin(), it);
					phi = out[index];
					psi = out[index + 1];
				}

				phi = numeric::conversions::degrees(phi);
				psi = numeric::conversions::degrees(psi);

				size_t index = bbdep_obj_sm.find_index_for_cdf_chi234(opt_vect[c].amino_acid, phi, psi);

				if(opt_vect[c].nchi == 1)
				{
					double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid, x[c], 0);
					out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));
				}
				else if(opt_vect[c].nchi == 2)
				{
					double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid, x[c], 0);
					out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

					if(first_line.r2 != 0)
					{
						double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            x[c + 1]);
						out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
					}
				}
				else if(opt_vect[c].nchi == 3)
				{
					double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid, x[c], 0);
					out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

					if(first_line.r2 != 0)
					{
						double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            x[c + 1]);
						out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
					}

					if(first_line.r3 != 0)
					{
						double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c + 1])),
						            x[c + 2]);
						out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
					}
				}
				else if(opt_vect[c].nchi == 4)
				{
					double c1 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chinumber(index, opt_vect[c].amino_acid, x[c], 0);
					out[c] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c1));

					if(first_line.r2 != 0)
					{
						double c2 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi1_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            x[c + 1]);
						out[c + 1] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c2));
					}

					if(first_line.r3 != 0)
					{
						double c3 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi12_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c + 1])),
						            x[c + 2]);
						out[c + 2] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c3));
					}

					if(first_line.r4 != 0)
					{
						double c4 = bbdep_obj_sm.get_degree_bbdep_from_phi_psi_x01_chi123_dep_actual_states(index,
						            opt_vect[c].amino_acid, numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c])),
						            numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c + 1])),
						            numeric::conversions::degrees(pepsgo::bbtools::to_positive_radians(out[c + 2])),
						            x[c + 3]);
						out[c + 3] = pepsgo::bbtools::normalize_to_mpi_to_ppi(numeric::conversions::radians(c4));
					}
				}
				c += opt_vect[c].nchi;
			}
		}
	}
}
