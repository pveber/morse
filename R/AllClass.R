# This file is part of morse
#
# morse is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# morse is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

setClass(Class = "survData",
         representation = representation(
           data_name = "character",
           data_df = "data.frame"))


setClass(Class = "survFit",
         representation = representation(
           model_name = "character", 
           model_pars = "character", 
           model_dims = "numeric",
           model_type = "character"))
