#This file is part of ChiVO, the Chilean Virtual Observatory
#A project sponsored by FONDEF (D11I1060)
#Copyright (C) 2015 Universidad Tecnica Federico Santa Maria Mauricio Solar
#                                                            Marcelo Mendoza
#                   Universidad de Chile                     Diego Mardones
#                   Pontificia Universidad Catolica          Karim Pichara
#                   Universidad de Concepcion                Ricardo Contreras
#                   Universidad de Santiago                  Victor Parada
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import numpy as np

def alinear(x1,y1,x2,y2,dim_1_x,dim_1_y,dim_2_x,dim_2_y, actual): # 1: Referencia 2: Pos. de la imagen actual
    diferencia_x = x1 - x2
    diferencia_y = y1 - y2

    matriz_final = np.zeros((dim_1_x,dim_1_y))

    for i in range(dim_2_x):
        for j in range(dim_2_y):
            x = i + diferencia_x
            y = j + diferencia_y

            if ( x > 0 and x < dim_1_x and y > 0 and y < dim_1_y):
                matriz_final[x][y] = actual[i][j]

    return matriz_final
