*/* RMNLIB - Library of useful routines for C and FORTRAN programming
* * Copyright (C) 1975-2012  Division de Recherche en Prevision Numerique
* *                          Environnement Canada
* *
* * This library is free software; you can redistribute it and/or
* * modify it under the terms of the GNU Lesser General Public
* * License as published by the Free Software Foundation,
* * version 2.1 of the License.
* *
* * This library is distributed in the hope that it will be useful,
* * but WITHOUT ANY WARRANTY; without even the implied warranty of
* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* * Lesser General Public License for more details.
* *
* * You should have received a copy of the GNU Lesser General Public
* * License along with this library; if not, write to the
* * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
* * Boston, MA 02111-1307, USA.
* */
*	FORTRAN to C bridging routines
*
*	RPN_COMM_chdir : bridge to f_RPN_COMM_chdir, null terminates string before calling C routine
*
        integer function RPN_COMM_chdir(string)
        character *(*)string
        integer f_RPN_COMM_chdir
        RPN_COMM_chdir=f_RPN_COMM_chdir(trim(string)//achar(0))
        return
        end
*
*       RPN_COMM_env_var : get value of environment variable
*       returns blanks if variable does not exist
*
	subroutine RPN_COMM_env_var(varname,value)
	implicit none
	character (len=*), intent(IN) :: varname
	character (len=*), intent(OUT) :: value

	integer status,length
        value = " "
        call RPN_COMM_getenv(trim(varname)//achar(0),value,len(value))
	return
	end
