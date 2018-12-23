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
#ifndef INCLUDED_get_dof_hh
#define INCLUDED_get_dof_hh

#include <core/pose/Pose.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/TorsionID.hh>

#include <vector>
#include <string>

#include "transform.hh"

namespace pepsgo
{

struct opt_element
{
    core::id::DOF_ID dofid;
    std::pair<core::Real, core::Real> bounds;
    core::chemical::AA amino_acid;
    std::string torsion_name;
    core::Size chainid;
    core::Size seqpos;
    core::Size nchi;
    core::Size ichi;

    opt_element();

    opt_element(core::id::DOF_ID _dofid, std::pair<core::Real, core::Real> _bounds, core::chemical::AA _amino_acid, std::string _torsion_name, core::Size _chainid,
                core::Size _seqpos,
                core::Size _nchi,
                core::Size _ichi);

    friend std::ostream& operator<<(std::ostream& stream, const opt_element& obj)
    {
        stream << obj.chainid << ' ' << obj.amino_acid << ' ' << obj.torsion_name << ' ' << obj.seqpos << ' ' << obj.nchi << ' ' << obj.ichi << ' ' << obj.bounds.first << ' ' << obj.bounds.second;
        return stream;
    }
};

std::vector<core::id::DOF_ID> get_peptide_sc_dof_THETA(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_sc_dof_THETA_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_sc_dof_D(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_sc_dof_D_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_sc_dof_PHI_without_CHI(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_sc_dof_PHI_without_CHI_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_THETA_without_mc(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_bb_dof_THETA_without_mc_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_D_without_mc(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_bb_dof_D_without_mc_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_without_mc_without_term(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);
std::vector<opt_element> get_peptide_bb_dof_PHI_without_mc_without_term_with_info(core::pose::Pose& pose, size_t start_ind, size_t stop_ind, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_first_term(core::pose::Pose& pose, int atom_type);
std::vector<opt_element> get_peptide_bb_dof_PHI_first_term_with_info(core::pose::Pose& pose, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_last_term(core::pose::Pose& pose, int atom_type);
std::vector<opt_element> get_peptide_bb_dof_PHI_last_term_with_info(core::pose::Pose& pose, int atom_type);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_term_chi1(core::pose::Pose& pose);
std::vector<opt_element> get_peptide_bb_dof_PHI_term_chi1_with_info(core::pose::Pose& pose);

std::vector<core::id::DOF_ID> get_peptide_bb_dof_PHI_term_chi2(core::pose::Pose& pose);
std::vector<opt_element> get_peptide_bb_dof_PHI_term_chi2_with_info(core::pose::Pose& pose);

std::vector<core::id::DOF_ID> get_peptide_mc_dof_THETA(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<opt_element> get_peptide_mc_dof_THETA_with_info(core::pose::Pose& pose, core::Size start, core::Size stop);

std::vector<core::id::DOF_ID> get_peptide_mc_dof_D(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<opt_element> get_peptide_mc_dof_D_with_info(core::pose::Pose& pose, core::Size start, core::Size stop);

std::vector<core::id::DOF_ID> get_all_dof_THETA(core::pose::Pose& pose);

std::vector<core::id::DOF_ID> get_all_dof_D(core::pose::Pose& pose);

std::vector<core::id::DOF_ID> get_phi_psi(core::pose::Pose& pose);
std::vector<opt_element> get_phi_psi_with_info(core::pose::Pose& pose);
std::vector<opt_element> get_phi_psi_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool remove_last);

std::vector<core::id::DOF_ID> get_omega(core::pose::Pose& pose);
std::vector<opt_element> get_omega_with_info(core::pose::Pose& pose);
std::vector<opt_element> get_omega_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool remove_last);

std::vector<core::id::DOF_ID> get_chi(core::pose::Pose& pose, core::Size start, core::Size stop);

std::vector<core::id::DOF_ID> get_chi1(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<core::id::DOF_ID> get_chi2(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<core::id::DOF_ID> get_chi3(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<core::id::DOF_ID> get_chi4(core::pose::Pose& pose, core::Size start, core::Size stop);

std::vector<core::id::DOF_ID> get_proline_chi1(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<opt_element> get_proline_chi1_with_info(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<core::id::DOF_ID> get_proline_chi2(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<opt_element> get_proline_chi2_with_info(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<core::id::DOF_ID> get_proline_chi3(core::pose::Pose& pose, core::Size start, core::Size stop);
std::vector<opt_element> get_proline_chi3_with_info(core::pose::Pose& pose, core::Size start, core::Size stop);

std::vector<std::tuple<core::id::DOF_ID, core::chemical::AA, core::Size, core::Size, core::Size> > get_peptide_all_chi_dof_with_aatype(core::pose::Pose& pose,
        core::Size start,
        core::Size stop);

std::vector<opt_element> get_peptide_all_chi_dof_with_info(core::pose::Pose& pose, core::Size start, core::Size stop, bool include_proline);
std::vector<opt_element> get_peptide_all_chi_dof_with_info(core::pose::Pose& pose, std::vector<core::Size> &ind, size_t chain_id, bool include_proline);


void insert_to_opt_vector(std::vector<opt_element> &opt_vector,
                          core::pose::Pose &pose,
                          std::string arguments,
                          transform::ranges &range);

}

#endif
