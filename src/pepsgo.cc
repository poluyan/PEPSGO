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
#include "pepsgo.hh"

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

namespace pepsgo
{
    
PEPSGO::PEPSGO()
{
    //score_fn = core::scoring::ScoreFunctionFactory::create_score_function("ref2015");
    score_fn = core::scoring::get_score_function();
    std::cout << "score function: " << score_fn->get_name() << std::endl;
}
void PEPSGO::set_peptide(std::string _peptide_sequence)
{
    peptide_sequence = _peptide_sequence;

    core::pose::make_pose_from_sequence(peptide, peptide_sequence,
                                        *core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));

    for(core::Size i = 1; i <= peptide.total_residue(); i++)
    {
        peptide.set_phi(i, -135.0);
        peptide.set_psi(i, 135.0);
        peptide.set_omega(i, 179.8);

        if(peptide.residue(i).type().aa() == core::chemical::aa_pro)
            continue;

        for(core::Size j = 1; j <= peptide.residue(i).nchi(); j++)
        {
            peptide.set_chi(j, i, 90.0);
        }
    }
    peptide.dump_pdb("output_pdb/peptide.pdb");
    
    ideal_peptide = peptide;

}
void PEPSGO::set_bbdep(std::string _bbdep_path)
{
    bbdep_path = _bbdep_path;
}
void PEPSGO::set_bbind(std::string _bbind_path)
{
    bbind_path = _bbind_path;
}

}
