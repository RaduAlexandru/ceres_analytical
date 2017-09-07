from sympy import *

init_printing(use_unicode=True)

#3d point
p3d_x, p3d_y, p3d_z=symbols('p3d_x p3d_y p3d_z')
p3d = Matrix([p3d_x, p3d_y, p3d_z])

#camera rotation, angle axis
w_x, w_y, w_z = symbols('w_x w_y w_z')
angle_axis=Matrix([w_x, w_y, w_z])
#camera translation
t_x, t_y, t_z = symbols('t_x t_y t_z')
T=Matrix([t_x, t_y, t_z])
#camera focal lengh
f = symbols('f')

#obs
obs_x, obs_y= symbols('obs_x obx_y')


#-----------------------------



#rotate
theta2 = angle_axis.dot(angle_axis);
theta = sqrt(theta2);
costheta = cos(theta);
sintheta = sin(theta);
theta_inverse = 1.0 / theta;
w=Matrix([w_x*theta_inverse, w_y*theta_inverse, w_z*theta_inverse])
w_cross_pt=w.cross(p3d)
tmp =(w[0] * p3d[0] + w[1] * p3d[1] + w[2] * p3d[2]) * ((1.0) - costheta);
p3d_rot=Matrix([p3d[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp,
              p3d[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp,
              p3d[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp])


#translate
p3d_rot_trans=Matrix([p3d_rot[0] + T[0],
                      p3d_rot[1] + T[1],
                      p3d_rot[2] + T[2]])

#dehomogenenous
p2d = Matrix([ -p3d_rot_trans[0] / p3d_rot_trans[2],
               -p3d_rot_trans[1] / p3d_rot_trans[2]])

#distorsion
#TODO
distorsion=1.0

#Compute final projected point position
predicted_x= f* p2d[0]
predicted_y= f* p2d[1]

#error
residual = Matrix([predicted_x - obs_x,
                   predicted_y - obs_y])



#jacobian with respect to rotation
j_w_x = residual.diff(w_x)
j_w_y = residual.diff(w_y)
j_w_z = residual.diff(w_z)
#jacbian translation
j_t_x = residual.diff(t_x)
j_t_y = residual.diff(t_y)
j_t_z = residual.diff(t_z)
#jacobian focal length
j_f = residual.diff(f)
#jacobian point3d
j_p3d_x = residual.diff(p3d_x)
j_p3d_y = residual.diff(p3d_y)
j_p3d_z = residual.diff(p3d_z)


##Normal
print "Normal"
print ("j_w_x")
print ((j_w_x))
print ("j_w_y")
print ((j_w_y))
print ("j_w_z")
print ((j_w_z))

print ("j_t_x")
print ((j_t_x))
print ("j_t_y")
print ((j_t_y))
print ("j_t_z")
print ((j_t_z))

print ("j_f")
print ((j_f))

print ("j_p3d_x")
print ((j_p3d_x))
print ("j_p3d_y")
print ((j_p3d_y))
print ("j_p3d_z")
print ((j_p3d_z))


##Simplyfied
print "Simplified"
print ("j_w_x")
print (simplify(j_w_x))
print ("j_w_y")
print (simplify(j_w_y))
print ("j_w_z")
print (simplify(j_w_z))

print ("j_t_x")
print (simplify(j_t_x))
print ("j_t_y")
print (simplify(j_t_y))
print ("j_t_z")
print (simplify(j_t_z))

print ("j_f")
print (simplify(j_f))

print ("j_p3d_x")
print (simplify(j_p3d_x))
print ("j_p3d_y")
print (simplify(j_p3d_y))
print ("j_p3d_z")
print (simplify(j_p3d_z))





# u, v, r = symbols('u v r')
#
# q = Matrix([u, v, sqrt(r*r - u*u - v*v)])
#
# q_u = q.diff(u)
# q_v = q.diff(v)
#
# E = q_u.dot(q_u)
# F = q_u.dot(q_v)
# G = q_v.dot(q_v)
#
# M_1 = Matrix([[E, F], [F, G]])
#
# #n = simplify(q_u.cross(q_v)) / (q_u.cross(q_v)).norm()
# n = Matrix([u/r, v/r, sqrt(r*r - u*u - v*v)/r])
#
#
# q_uu = simplify(q_u.diff(u))
# q_vv = simplify(q_v.diff(v))
# q_uv = simplify(q_u.diff(v))
#
#
# L = simplify(q_uu.dot(n))
# M = simplify(q_uv.dot(n))
# N = simplify(q_vv.dot(n))
#
#
# M_2 = Matrix(2, 2, [L, M, M, N])
#
# # ===========================
#
# k_G = simplify((L * N - M**2) / (E * G - F**2))
# k_M = simplify((1/2) * (E * N + G * L - 2 * F * M) / (E * G - F**2))
# k = var('k')
# f = k**2 - 2 * k_M * k + k_G
# prin=solve(f, k)
# print (latex(prin))
# print("\kappa_G = " + latex(k_G) + ", 'kappa_M = " + latex(k_M) + ",\kappa = " + latex(solve(f, k)))
