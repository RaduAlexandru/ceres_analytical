from sympy import *

init_printing(use_unicode=True)

#3d point
p3d_x, p3d_y, p3d_z=symbols('p3d_x p3d_y p3d_z')
p3d = Matrix([p3d_x, p3d_y, p3d_z])

#camera rotation, angle axis
q_w, q_x, q_y, q_z = symbols('q_w q_x q_y q_z')
q=Matrix([q_w, q_x, q_y, q_z ])
#camera translation
t_x, t_y, t_z = symbols('t_x t_y t_z')
T=Matrix([t_x, t_y, t_z])
#camera focal lengh
f = symbols('f')

#obs
obs_x, obs_y= symbols('obs_x obx_y')


#-----------------------------

#get unit quat
scale = 1.0 / sqrt( q[0] * q[0] +
                    q[1] * q[1] +
                    q[2] * q[2] +
                    q[3] * q[3]);
u_q=Matrix([scale * q[0],
            scale * q[1],
            scale * q[2],
            scale * q[3]])



#rotate
t2 =  u_q[0] * u_q[1];
t3 =  u_q[0] * u_q[2];
t4 =  u_q[0] * u_q[3];
t5 = -u_q[1] * u_q[1];
t6 =  u_q[1] * u_q[2];
t7 =  u_q[1] * u_q[3];
t8 = -u_q[2] * u_q[2];
t9 =  u_q[2] * u_q[3];
t1 = -u_q[3] * u_q[3];
p3d_rot=Matrix([ 2* ((t8 + t1) * p3d[0] + (t6 - t4) * p3d[1] + (t3 + t7) * p3d[2]) + p3d[0],
                 2 * ((t4 + t6) * p3d[0] + (t5 + t1) * p3d[1] + (t9 - t2) * p3d[2]) + p3d[1],
                 2 * ((t7 - t3) * p3d[0] + (t2 + t9) * p3d[1] + (t5 + t8) * p3d[2]) + p3d[2]])



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
j_q_w = residual.diff(q_w)
j_q_x = residual.diff(q_x)
j_q_y = residual.diff(q_y)
j_q_z = residual.diff(q_z)
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
print ("j_q_w")
print ((j_q_w))
print ("j_q_x")
print ((j_q_x))
print ("j_q_y")
print ((j_q_y))
print ("j_q_z")
print ((j_q_z))

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
print ("j_q_w")
print (simplify(j_q_w))
print ("j_q_x")
print (simplify(j_q_x))
print ("j_q_y")
print (simplify(j_q_y))
print ("j_q_z")
print (simplify(j_q_z))

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
