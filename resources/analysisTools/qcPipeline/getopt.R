#!/usr/bin/env Rscript
#
# This file is part of the AlignmentAndQCWorkflow plugin.
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 or 3 of the License.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this script.  If not, see <https://www.gnu.org/licenses/>.
#

require(getopt)
source(file.path(dirname(get_Rscript_filename()), "qq.R"));

getopt2 = function(..., attach = parent.frame()) {
	global_obj_name = ls(envir = attach)
	
	m = list(...)[[1]]
	args = commandArgs(TRUE)
	
	if(length(args) == 0) {
		cat(getopt::getopt(m, usage = TRUE))
		q(save = "no")
	}
	
	for(i in seq_len(nrow(m))) {
		if(length(grep("^[a-zA-Z\\_.][0-9a-zA-Z\\_.]*$", m[i, 1])) == 0) {
			cat(qq("error: long option name [@{m[i, 1]}] is not allowed. Option long name can only be '^[a-zA-Z\\_.][0-9a-zA-Z\\_.]*$'\n\n\n"))
			cat(getopt::getopt(m, usage = TRUE))
			q(save = "no")
		}
	}
	
	for(i in seq_len(nrow(m))) {
		if(m[i, 3] == "1") {
			if(length(grep(qq("^-{1,2}@{m[i, 1]}$"), args)) == 0 && 
			   length(grep(qq("^-{1,2}@{m[i, 2]}$"), args)) == 0) {
				
				cat(qq("error: [--@{m[i, 1]}|-@{m[i, 2]}] is a mandatory field, but you have not specified it.\n\n\n"))
				cat(getopt::getopt(m, usage = TRUE))
				q(save = "no")
			}
		}
	}
	
	opt = getopt::getopt(...)

	opt_obj_name = ls(envir = as.environment(opt))
	
	# specified from command line
	specified_opt_obj_name = opt_obj_name[!is.null(opt_obj_name)]
	# export to global environment
	for(o in specified_opt_obj_name) {
		assign(o, opt[[o]], envir = attach)
	}

	# defined with default values while not specified in command line
	specified_global_obj_name = intersect(opt_obj_name[is.null(opt_obj_name)], global_obj_name)
	# already have, do nothing
	
	# undefined values
	rest_opt_obj_name = setdiff(opt_obj_name[is.null(opt_obj_name)], global_obj_name)
	for(o in rest_opt_obj_name) {
		warning(qq("@{o} has not be specified in command line and it does not have a default value either.\n"));
	}
	return(invisible(NULL))
	
}
