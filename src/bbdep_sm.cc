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
#include "bbdep_sm.hh"
#include "bbutils.hh"

namespace pepsgo
{
namespace bbdep
{

BBDEP_Dunbrack_sm::BBDEP_Dunbrack_sm(std::string _path_to_files, size_t _cdf_grid_step):
    path_to_files(_path_to_files), cdf_grid_step(_cdf_grid_step) {}

// amino acids in 1letter code format
void BBDEP_Dunbrack_sm::initialize_all(bool create_cdf_sum, std::string amino_acids)
{
    load_data_em();
}


void BBDEP_Dunbrack_sm::load_data_em()
{
    aa_sm_1d.resize(4); // 0 ser, 1 val, 2 cys, 3 thr
    pepsgo::bbdep::load_data_sm(path_to_files + "ser.bbdep.rotamers.lib.gz", aa_sm_1d[0].lib, aa_sm_1d[0].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "val.bbdep.rotamers.lib.gz", aa_sm_1d[1].lib, aa_sm_1d[1].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "cys.bbdep.rotamers.lib.gz", aa_sm_1d[2].lib, aa_sm_1d[2].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "thr.bbdep.rotamers.lib.gz", aa_sm_1d[3].lib, aa_sm_1d[3].libn);
    
    aa_sm_2d.resize(8); // 0 trp, 1 his, 2 asn, 3 asp, 4 phe, 5 tyr, 6 ile, 7 leu
    pepsgo::bbdep::load_data_sm(path_to_files + "trp.bbdep.rotamers.lib.gz", aa_sm_2d[0].lib, aa_sm_2d[0].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "his.bbdep.rotamers.lib.gz", aa_sm_2d[1].lib, aa_sm_2d[1].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "asn.bbdep.rotamers.lib.gz", aa_sm_2d[2].lib, aa_sm_2d[2].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "asp.bbdep.rotamers.lib.gz", aa_sm_2d[3].lib, aa_sm_2d[3].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "phe.bbdep.rotamers.lib.gz", aa_sm_2d[4].lib, aa_sm_2d[4].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "tyr.bbdep.rotamers.lib.gz", aa_sm_2d[5].lib, aa_sm_2d[5].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "ile.bbdep.rotamers.lib.gz", aa_sm_2d[6].lib, aa_sm_2d[6].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "leu.bbdep.rotamers.lib.gz", aa_sm_2d[7].lib, aa_sm_2d[7].libn);
    
    aa_sm_3d.resize(3); // 0 met, 1 glu, 2 gln, 3 pro
    pepsgo::bbdep::load_data_sm(path_to_files + "met.bbdep.rotamers.lib.gz", aa_sm_3d[0].lib, aa_sm_3d[0].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "glu.bbdep.rotamers.lib.gz", aa_sm_3d[1].lib, aa_sm_3d[1].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "pro.bbdep.rotamers.lib.gz", aa_sm_3d[2].lib, aa_sm_3d[2].libn);
    
    aa_sm_4d.resize(2); // 0 arg, 1 lys
    pepsgo::bbdep::load_data_sm(path_to_files + "arg.bbdep.rotamers.lib.gz", aa_sm_4d[0].lib, aa_sm_4d[0].libn);
    pepsgo::bbdep::load_data_sm(path_to_files + "lys.bbdep.rotamers.lib.gz", aa_sm_4d[1].lib, aa_sm_4d[1].libn);
}


}
}
