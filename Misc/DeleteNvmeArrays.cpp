/*
 *
 * Copyright 2018 The RMG Project Developers. See the COPYRIGHT file 
 * at the top-level directory of this distribution or in the current
 * directory.
 * 
 * This file is part of RMG. 
 * RMG is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * any later version.
 *
 * RMG is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "main.h"


// Cleans up any mmapped NVME arrays we may have created
void DeleteNvmeArrays(void)
{
    std::string newpath;

    if(ct.nvme_orbital_fd > 0)
    {
        ftruncate(ct.nvme_orbital_fd, 0);
        close(ct.nvme_orbital_fd);
        newpath = ct.nvme_orbitals_path + std::string("rmg_orbital") + std::to_string(pct.spinpe) + "_" +
                  std::to_string(pct.kstart) + "_" + std::to_string(pct.gridpe);
        unlink(newpath.c_str());
    }

    if(ct.nvme_work_fd > 0)
    {
        ftruncate(ct.nvme_work_fd, 0);
        close(ct.nvme_work_fd);
        newpath = ct.nvme_work_path + std::string("rmg_work") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        unlink(newpath.c_str());
    }

    if(ct.nvme_weight_fd > 0)
    {
        ftruncate(ct.nvme_weight_fd, 0);
        close(ct.nvme_weight_fd);
        newpath = ct.nvme_weights_path + std::string("rmg_weight") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        unlink(newpath.c_str());
    }

    if(ct.nvme_Bweight_fd > 0)
    {
        ftruncate(ct.nvme_Bweight_fd, 0);
        close(ct.nvme_Bweight_fd);
        newpath = ct.nvme_weights_path + std::string("rmg_Bweight") + std::to_string(pct.spinpe) +
                  std::to_string(pct.kstart) + std::to_string(pct.gridpe);
        unlink(newpath.c_str());
    }

}
