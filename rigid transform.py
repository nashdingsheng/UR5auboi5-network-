# -*- coding: utf-8 -*-
"""
Created on Mar 2018

@author: March
"""
# import xmlrpclib
# from SimpleXMLRPCServer import SimpleXMLRPCServer
# import urlib

import math


def calculate(pt1, pt2, pt3):
    np = numpy()
    points = [pt1[0:3], pt2[0:3], pt3[0:3]]
    origin, point1, point2 = np.sub(points, points[0])
    para_plane = np.dot(np.inv(points), [[-1], [-1], [-1]])
    para_plane = [para_plane[0][0], para_plane[1][0], para_plane[2][0]]
    axis_y = np.div(point1, np.norm(point1))
    x_left = [axis_y[0:2], para_plane[0:2]]
    if para_plane[0:2] == [0, 0]:
        c1 = -1 / para_plane[2] - points[0][2]          # a1, b1, c1 is normal vector of axis X
        a0, b0, c0 = axis_y                             # a0, b0, c0 is normal vector of axis Y
        a = a0 ** 2 + b0 ** 2
        b = 2 * b0 * c0 * c1                            # a*x**2 + b*x +c =0,  x = b0
        c = a0 * a0 * c1 * c1 - a0 ** 2 + c0 * c0 * c1 * c1
        b10 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        b11 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
        if a0 == 0:
            a10 = 1
            a11 = -1
        else:
            a10 = (-b0 * b10 - c0 * c1) / a0
            a11 = (-b0 * b11 - c0 * c1) / a0
        if a10 ** 2 + b10 ** 2 + c1 ** 2 == 1:
            axis_x = [a10, b10, c1]
        else:
            axis_x = [a11, b11, c1]
    elif x_left[0][0] == 0.0 and x_left[1][0] == 0.0:
        a0, b0, c0 = axis_y
        c1 = (1 + para_plane[1] * points[0][1] + para_plane[2] * points[0][2])/(c0 * para_plane[1] / b0 - para_plane[2])
        b1 = -(c0 * c1 / b0)
        axis_x = np.div([1, b1, c1], np.norm([1, b1, c1]))
    else:
        axis_x_ab = np.dot(np.inv(x_left), [[-axis_y[2]],
                                            [(-1) * para_plane[2] - 1 - np.dot(para_plane, points[0])]])
        axis_x = np.div([axis_x_ab[0][0], axis_x_ab[1][0], 1], np.norm([axis_x_ab[0][0], axis_x_ab[1][0], 1]))
    axis_z = np.cross(axis_x, axis_y)
    pos_z = np.inv(np.rot_vec2mat(pt1[3:6]))[2]
    print('*************************', pos_z, axis_z, np.dot(axis_z, pos_z))

    if np.dot(axis_z, pos_z) < 0:
        axis_x = np.mul(-1, axis_x)
        axis_z = np.cross(axis_x, axis_y)
    ang = np.cal_ang(axis_x, axis_y, axis_z)


    mt =[axis_x, axis_y, axis_z]
    mtt = [[axis_x[0],axis_y[0],axis_z[0]], [axis_x[1],axis_y[1],axis_z[1]], [axis_x[2],axis_y[2],axis_z[2]]]
    rx, ry, rz = np.rot_vector(mtt)
    # mt2 = np.dot(np.rotate_eta([np.deg2rad(90),np.deg2rad(0),np.deg2rad(0)]), mt)
    tx, ty, tz = np.rot_vector(mt)

    axis_x.extend(axis_y)
    axis_x.extend(axis_z)
    return axis_x, [rx, ry, rz, tx, ty, tz], ang


def transform(mat, pnt):
    np = numpy()
    axis_x = mat[0:3]
    axis_y = mat[3:6]
    axis_z = mat[6:9]
    return [np.dot(pnt, axis_x), np.dot(pnt, axis_y), np.dot(pnt, axis_z)]


def transform_inv(mat, pnt):
    np = numpy()
    mt = [mat[0:3], mat[3:6], mat[6:9]]
    axis_x, axis_y, axis_z = np.inv(mt)

    return [np.dot(pnt, axis_x), np.dot(pnt, axis_y), np.dot(pnt, axis_z)]


class rotation(object):
    def __init__(self):
        self.axis_x = [1, 0, 0]
        self.axis_y = [0, 1, 0]
        self.axis_z = [0, 0, 1]
        self.theta = [0, 0, 0]
        self.eta = [0, 0, 0]

    def calculation(self, points):
        np = numpy()
        origin, point1, point2 = np.sub(points, points[0])
        para_plane = np.dot(np.inv(points), [[-1], [-1], [-1]])
        para_plane = [para_plane[0][0], para_plane[1][0], para_plane[2][0]]
        self.axis_y = np.div(point1, np.norm(point1))
        x_left = [self.axis_y[0:2], para_plane[0:2]]
        if para_plane[0:2] == [0, 0]:
            c1 = -1 / para_plane[2] - points[0][2]     # a1, b1, c1 is normal vector of axis Y
            a0, b0, c0 = self.axis_y                   # a0, b0, c0 is normal vector of axis X
            a = a0**2 + b0**2
            b = 2 * b0 * c0 * c1                       # a*x**2 + b*x +c =0,  x = b1
            c = a0 * a0 * c1 * c1 - a0**2 + c0 * c0 * c1 * c1
            b10 = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
            b11 = (-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            a10 = (-b0 * b10 - c0 * c1) / a0
            a11 = (-b0 * b11 - c0 * c1) / a0
            if a10**2 + b10**2 + c1**2 == 1:
                self.axis_x = [a10, b10, c1]
            else:
                self.axis_x = [a11, b11, c1]
        else:
            axis_x_ab = np.dot(np.inv(x_left), [[-self.axis_y[2]],
                                                      [(-1) * para_plane[2] - 1 - np.dot(para_plane, points[0])]])
            self.axis_x = np.div([axis_x_ab[0][0], axis_x_ab[1][0], 1], np.norm([axis_x_ab[0][0], axis_x_ab[1][0], 1]))
        self.axis_z = np.cross(self.axis_x, self.axis_y)
        if np.dot(self.axis_z, [0, 0, 1]) < 0:
            self.axis_x = np.mul(-1, self.axis_x)
            self.axis_z = np.cross(self.axis_x, self.axis_y)

        # if self.axis_z[0] > 0 and self.axis_z[1] < 0 and self.axis_z[2] < 0:
        #     self.axis_x = np.mul(self.axis_x, -1)
        #     self.axis_z = np.cross(self.axis_x, self.axis_y)

        theta_x = math.atan(-self.axis_y[2] / self.axis_z[2])
        theta_y = math.asin(self.axis_x[2])
        theta_z = math.atan(-self.axis_x[1] / self.axis_x[0])
        if self.axis_x[2] >= 0 and self.axis_z[2] < 0:
            theta_x += math.pi
        if self.axis_x[2] >= 0 and self.axis_x[0] < 0:
            theta_z += math.pi
        self.theta = [theta_x, theta_y, theta_z]

        eta_y = math.asin(-self.axis_z[0])
        eta_x = math.atan(self.axis_z[1] / self.axis_z[2])  # np.arcsin(axis_z[1]/np.cos(eta_y))
        eta_z = math.atan(self.axis_y[0] / self.axis_x[0])  # np.arcsin(axis_y[0]/np.cos(eta_y))
        if self.axis_x[0] < 0 and self.axis_y[0] > 0:
            eta_z += math.pi
        self.eta = [eta_x, eta_y, eta_z]

        return self.theta

    def get_eta(self):
        return self.eta

    def transform(self, pnt):
        np = numpy()
        return [np.dot(pnt, self.axis_x), np.dot(pnt, self.axis_y), np.dot(pnt, self.axis_z)]


class numpy(object):
    def add(self, list1, list2):
        if type(list1) == list and type(list2) == list:
            if type(list1[0]) == list and type(list2[0]) == list:
                result = []
                for a, b in zip(list1, list2):
                    result.append([c + d for c, d in zip(a, b)])
                return result
            elif type(list1[0]) == list and type(list2[0]) != list:
                result = []
                for a in list1:
                    result.append([c + d for c, d in zip(a, list2)])
                return result
            elif type(list1[0]) != list and type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c + d for c, d in zip(a, list1)])
                return result
            else:
                return [a + b for a, b in zip(list1, list2)]
        elif type(list1) == list and type(list2) != list:
            if type(list1[0]) == list:
                result = []
                for a in list1:
                    result.append([c + list2 for c in a])
                return result
            else:
                return [a + list2 for a in list1]
        elif type(list1) != list and type(list2) == list:
            if type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c + list1 for c in a])
                return result
            else:
                return [a + list1 for a in list2]

    def sub(self, list1, list2):
        if type(list1) == list and type(list2) == list:
            if type(list1[0]) == list and type(list2[0]) == list:
                result = []
                for a, b in zip(list1, list2):
                    result.append([c - d for c, d in zip(a, b)])
                return result
            elif type(list1[0]) == list and type(list2[0]) != list:
                result = []
                for a in list1:
                    result.append([c - d for c, d in zip(a, list2)])
                return result
            elif type(list1[0]) != list and type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c - d for c, d in zip(a, list1)])
                return result
            else:
                return [a - b for a, b in zip(list1, list2)]
        elif type(list1) == list and type(list2) != list:
            if type(list1[0]) == list:
                result = []
                for a in list1:
                    result.append([c - list2 for c in a])
                return result
            else:
                return [a - list2 for a in list1]
        elif type(list1) != list and type(list2) == list:
            if type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c - list1 for c in a])
                return result
            else:
                return [a - list1 for a in list2]

    def mul(self, list1, list2):
        if type(list1) == list and type(list2) == list:
            if type(list1[0]) == list and type(list2[0]) == list:
                result = []
                for a, b in zip(list1, list2):
                    result.append([c * d for c, d in zip(a, b)])
                return result
            elif type(list1[0]) == list and type(list2[0]) != list:
                result = []
                for a in list1:
                    result.append([c * d for c, d in zip(a, list2)])
                return result
            elif type(list1[0]) != list and type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c * d for c, d in zip(a, list1)])
                return result
            else:
                return [a * b for a, b in zip(list1, list2)]
        elif type(list1) == list and type(list2) != list:
            if type(list1[0]) == list:
                result = []
                for a in list1:
                    result.append([c * list2 for c in a])
                return result
            else:
                return [a * list2 for a in list1]
        elif type(list1) != list and type(list2) == list:
            if type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c * list1 for c in a])
                return result
            else:
                return [a * list1 for a in list2]

    def div(self, list1, list2):
        if type(list1) == list and type(list2) == list:
            if type(list1[0]) == list and type(list2[0]) == list:
                result = []
                for a, b in zip(list1, list2):
                    result.append([c / d for c, d in zip(a, b)])
                return result
            elif type(list1[0]) == list and type(list2[0]) != list:
                result = []
                for a in list1:
                    result.append([c / d for c, d in zip(a, list2)])
                return result
            elif type(list1[0]) != list and type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c / d for c, d in zip(a, list1)])
                return result
            else:
                return [a / b for a, b in zip(list1, list2)]
        elif type(list1) == list and type(list2) != list:
            if type(list1[0]) == list:
                result = []
                for a in list1:
                    result.append([c / list2 for c in a])
                return result
            else:
                return [a / list2 for a in list1]
        elif type(list1) != list and type(list2) == list:
            if type(list2[0]) == list:
                result = []
                for a in list2:
                    result.append([c / list1 for c in a])
                return result
            else:
                return [a / list1 for a in list2]

    def rad2deg(self, rads):
        if type(rads) != list:
            return rads * 180 / math.pi
        else:
            return [rad * 180 / math.pi for rad in rads]

    def deg2rad(self, angls):
        if type(angls) != list:
            return angls * math.pi / 180
        else:
            return [ang * math.pi / 180 for ang in angls]

    def rotate_eta(self, angles):
        (c1, s1), (c2, s2), (c3, s3) = [(math.cos(a), math.sin(a)) for a in angles]
        r_x = [[1, 0, 0], [0, c1, -s1], [0, s1, c1]]
        r_y = [[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]]
        r_z = [[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]]

        return self.dot(r_z, self.dot(r_y, r_x))

    def rotate_theta(self, angles):
        (c1, s1), (c2, s2), (c3, s3) = [(math.cos(a), math.sin(a)) for a in angles]
        r_x = [[1, 0, 0], [0, c1, -s1], [0, s1, c1]]
        r_y = [[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]]
        r_z = [[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]]

        return self.dot(r_x, self.dot(r_y, r_z))

    def cross(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0]]

        return c

    def dot(self, pa, pb):
        if type(pa[0]) == list:
            rows_a = len(pa)
            cols_a = len(pa[0])
            if type(pb[0]) == list:
                rows_b = len(pb)
                cols_b = len(pb[0])
            else:
                rows_b = len(pb)
                cols_b = 1

            if cols_a != rows_b:
              print("Cannot multiply the two matrices. Incorrect dimensions.")
              return
            # Create the result matrix
            # Dimensions would be rows_A x cols_B
            C = [[0 for row in range(cols_b)] for col in range(rows_a)]
            for i in range(rows_a):
                for j in range(cols_b):
                    for k in range(cols_a):
                        C[i][j] += pa[i][k] * pb[k][j]
            return C
        else:
            result = 0
            for a, b in zip(pa, pb):
                result += a * b
            return result

    def norm(self, point):
        x, y, z = point
        return math.sqrt(x**2 + y**2 + z**2)

    def intp(self, arr):
        return [int(a*10) / 10 for a in arr]

    def inv(self, A):
        if len(A) == 3:
            detA = A[0][0]*A[1][1]*A[2][2] + A[1][0]*A[2][1]*A[0][2] + A[2][0]*A[0][1]*A[1][2] \
                   - A[0][0]*A[2][1]*A[1][2] - A[2][0]*A[1][1]*A[0][2] - A[1][0]*A[0][1]*A[2][2]
            a00 = A[1][1] * A[2][2] - A[1][2] * A[2][1]
            a01 = A[0][2] * A[2][1] - A[0][1] * A[2][2]
            a02 = A[0][1] * A[1][2] - A[0][2] * A[1][1]
            a10 = A[1][2] * A[2][0] - A[1][0] * A[2][2]
            a11 = A[0][0] * A[2][2] - A[0][2] * A[2][0]
            a12 = A[0][2] * A[1][0] - A[0][0] * A[1][2]
            a20 = A[1][0] * A[2][1] - A[1][1] * A[2][0]
            a21 = A[0][1] * A[2][0] - A[0][0] * A[2][1]
            a22 = A[0][0] * A[1][1] - A[0][1] * A[1][0]
            if detA !=0:
                result = [[a00 / detA, a01 / detA, a02 / detA],
                          [a10 / detA, a11 / detA, a12 / detA],
                          [a20 / detA, a21 / detA, a22 / detA]]

            return result
        else:
            detA = A[0][0] * A[1][1] - A[0][1] * A[1][0]
            return [[A[1][1] / detA, -A[0][1] / detA], [-A[1][0] / detA, A[0][0] / detA]]

    def rot_vector(self, mat):
        r = [mat[2][1] - mat[1][2], mat[0][2] - mat[2][0], mat[1][0] - mat[0][1]]
        s = math.sqrt((r[0] ** 2 + r[1] ** 2 + r[2] ** 2) * 0.25)
        c = (mat[0][0] + mat[1][1] + mat[2][2] - 1) * 0.5
        c = min(max(c, -1), 1)
        theta = math.acos(c)
        if s < 1e-5:
            if c > 0:
                r = [0, 0, 0]
            else:
                t = (mat[0][0] + 1) * 0.5
                r[0] = math.sqrt(max(t, 0.0))
                t = (mat[1][1] + 1) * 0.5
                r[1] = math.sqrt(max(t, 0.0)) * (-1 if mat[0][1] < 0 else 1)
                t = (mat[2][2] + 1) * 0.5
                r[2] = math.sqrt(max(t, 0.0)) * (-1 if mat[0][2] < 0 else 1)
                if math.fabs(r[0]) < math.fabs(r[1]) and math.fabs(r[0]) < math.fabs(r[2]) \
                        and (mat[1][2] > 0) != (r[1] * r[2] > 0):
                    r[2] = -r[2]
                theta /= math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                r = self.mul(r, theta)
        else:
            vth = 1 / (2 * s)
            vth *= theta
            r = self.mul(r, vth)

        return r

    def rot_vec2mat(self, r):
        ang = self.norm(r)
        rx, ry, rz = self.div(r, ang)
        a = self.mul(math.cos(ang),[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        b = self.mul(1 - math.cos(ang),
                     [[rx * rx, rx * ry, rx * rz], [ry * rx, ry * ry, ry * rz], [rz * rx, rz * ry, rz * rz]])
        c = self.mul(math.sin(ang), [[0, -rz, ry], [rz, 0, -rx], [-ry, rx, 0]])

        return self.add(a, self.add(b, c))

    def cal_ang(self, ax, ay, az):
        def check(xx, yy, zz):
            check_ax, check_ay, check_az = self.rotate_theta([xx, yy, zz])
            return self.intp(self.sub(check_ax, ax)) == [0.0, 0.0, 0.0] and \
                   self.intp(self.sub(check_ay, ay)) == [0.0, 0.0, 0.0] and \
                   self.intp(self.sub(check_az, az)) == [0.0, 0.0, 0.0]

        ry = math.asin(ax[2])
        if ax[0] == 0.0 and az[2] != 0.0:
            if 0 < -ax[1] / math.cos(ry):
                tempt = -ax[1] / math.cos(ry) if -ax[1] / math.cos(ry) < 1.0 else 1.0
            else:
                tempt = -ax[1] / math.cos(ry) if -ax[1] / math.cos(ry) > -1.0 else -1.0
            rz = math.asin(tempt)
            rx = math.atan(-ay[2] / az[2])
        elif az[2] == 0.0 and ax[0] != 0.0:
            rx = math.asin(-ay[2] / math.cos(ry))
            rz = math.atan(-ax[1] / ax[0])
        elif ax[0] == 0.0 and az[2] == 0.0:
            rz = math.asin(-ax[1] / math.cos(ry))
            rx = math.asin(-ay[2] / math.cos(ry))
        else:
            rx = math.atan(-ay[2] / az[2])
            rz = math.atan(-ax[1] / ax[0])

        if ax[2] >= 0:
            if ax[0] >= 0 and az[2] >= 0:
                if check(rx, ry, rz):
                    rx += 0
            elif ax[0] >= 0 and az[2] < 0:
                if check(rx + math.pi, ry, rz) and check(rx, math.pi - ry, rz + math.pi):
                    if math.fabs(rx + math.pi) + math.fabs(ry) + math.fabs(rz) > \
                                            math.fabs(rx) + math.fabs(math.pi - ry) + math.fabs(rz + math.pi):
                        ry = math.pi - ry
                        rz += math.pi
                    else:
                        rx += math.pi
            elif ax[0] < 0 and az[2] >= 0:
                if check(rx, ry, rz + math.pi) and check(rx + math.pi, math.pi - ry, rz):
                    if math.fabs(rx) + math.fabs(ry) + math.fabs(rz + math.pi) > \
                                            math.fabs(rx + math.pi) + math.fabs(math.pi - ry) + math.fabs(rz):
                        ry = math.pi - ry
                        rx += math.pi
                    else:
                        rz += math.pi
            elif ax[0] < 0 and az[2] < 0:
                if check(rx + math.pi, ry, rz + math.pi) and check(rx, math.pi - ry, rz):
                    if math.fabs(rx + math.pi) + math.fabs(ry) + math.fabs(rz + math.pi) > \
                                            math.fabs(rx) + math.fabs(math.pi - ry) + math.fabs(rz):
                        ry = math.pi - ry
                    else:
                        rx += math.pi
                        rz += math.pi
                elif check(rx + math.pi, ry, rz + math.pi):
                    rx += math.pi
                    rz += math.pi
                elif check(rx, math.pi - ry, rz):
                    ry = math.pi - ry
        else:
            if ax[0] >= 0 and az[2] >= 0:
                if check(rx, ry, rz):
                    rx += 0
            elif ax[0] >= 0 and az[2] < 0:
                if check(rx + math.pi, ry, rz) and check(rx, math.pi - ry, rz + math.pi):
                    if math.fabs(rx + math.pi) + math.fabs(ry) + math.fabs(rz) > \
                                            math.fabs(rx) + math.fabs(math.pi - ry) + math.fabs(rz + math.pi):
                        ry = math.pi - ry
                        rz += math.pi
                    else:
                        rx += math.pi
            elif ax[0] < 0 and az[2] >= 0:
                if check(rx, ry, rz + math.pi) and check(rx + math.pi, math.pi - ry, rz):
                    if math.fabs(rx) + math.fabs(ry) + math.fabs(rz + math.pi) > \
                                            math.fabs(rx + math.pi) + math.fabs(math.pi - ry) + math.fabs(rz):
                        ry = math.pi - ry
                        rx += math.pi
                    else:
                        rz += math.pi
            elif ax[0] < 0 and az[2] < 0:
                if check(rx + math.pi, ry, rz + math.pi) and check(rx, -math.pi - ry, rz):
                    if math.fabs(rx + math.pi) + math.fabs(ry) + math.fabs(rz + math.pi) > \
                                            math.fabs(rx) + math.fabs(-math.pi - ry) + math.fabs(rz):
                        ry = -math.pi - ry
                    else:
                        rx += math.pi
                        rz += math.pi
                if check(rx + math.pi, ry, rz + math.pi):
                    rx += math.pi
                    rz += math.pi
                elif check(rx, -math.pi - ry, rz):
                    ry = -math.pi - ry
        if rx > math.pi:
            rx = rx - math.pi * 2
        if ry > math.pi:
            ry = ry - math.pi * 2
        if rz > math.pi:
            rz = rz - math.pi * 2
        if rx < -math.pi:
            rx = rx + math.pi * 2
        if ry < -math.pi:
            ry = ry + math.pi * 2
        if rz < -math.pi:
            rz = rz + math.pi * 2

        return [rx, ry, rz]

def rpy2rotate(rpy):
    rx,ry,rz =np.array(rpy)
    matx = [[1,0,0],[0,np.cos(rx),np.sin(rx)],[0,-np.sin(rx), np.cos(rx)]]
    maty = [[np.cos(ry), 0, -np.sin(ry)],[0,1,0],[np.sin(ry),0,np.cos(ry)]]
    matz = [[np.cos(rz), np.sin(rz),0],[-np.sin(rz), np.cos(rz),0],[0,0,1]]
    tranmartix = np.dot(matx,np.dot(maty,matz))#这里需要确认欧拉角的旋转顺序是Z,Y,X轴
    tranmartix =np.linalg.inv(tranmartix)
    line_unit_vector = tranmartix[2:3,:][0]
    return tranmartix
def homo2quater(homo):
    ortho_base = homo.flatten()
    w = np.sqrt(1 + ortho_base[0] + ortho_base[4] + ortho_base[8]) / 2;
    x = (ortho_base[7] - ortho_base[5]) / (4 * w);
    y = (ortho_base[2] - ortho_base[6]) / (4 * w);
    z = (ortho_base[3] - ortho_base[1]) / (4 * w);
    return (w,x,y,z)
def quater2homo(quater_para):
    w,x,y,z = quater_para[0],quater_para[1],quater_para[2],quater_para[3]
    xx = x*x
    yy = y*y
    zz = z*z
    xy = x*y
    wz = w*z
    wy = w*y
    xz = x*z
    yz = y*z
    wx = w*x
    homo = [[1.0-2*(yy+zz),2*(xy-wz),2*(wy+xz)],
            [2*(xy+wz),1.0-2*(xx+zz),2*(yz-wx)],
            [2*(xz-wy),2*(yz+wx),1.0-2*(xx+yy)]]  
    return homo







C++  旋转矩阵 Eigen::Matrix 对应的Matrix3d/3f等等 ，四元数 Eigen::Quaternion<_Scalar, _Options>，旋转矢量 Eigen::Angle_Axis<_Scaler>，
      之间的转换可以参考各自类对应的构造函数以及成员函数就能完成所以的转换。
    由于欧拉角的顺序问题，在Eigen::Geometry中没有很多介绍欧拉角转换，因为当加入顺序自由度会很麻烦，所以留给我们自己进行确立：
     关于欧拉角转的一种情况(Eigen::Geometry)给出一个简单例子,以后旋转轴顺序可以自己定义：
     
    Matrix3f m;
    m = AngleAxisf(angle1, Vector3f::UnitZ())
    AngleAxisf(angle2, Vector3f::UnitY())
    AngleAxisf(angle3, Vector3f::UnitZ());
     Eigen::Geometry中，还出了一种AngleaAxis类， 其实也是比较好的。 很清晰挣了欧拉角问题，最主要的是在实际的空间几何中，能更加形象表示Affine变换。

          
关于平面上一个向量围绕着某一个固定点进行旋转的问题，图片参考 robot_model.jpg(everything search)        
        
def CB_path_planning(p1,p2,start_point,theta_rotate):
    direction_vector = np.array(p2)-np.array(p1)
    direction_vector = direction_vector/np.linalg.norm(direction_vector)
    quaternion_para =[np.cos(theta_rotate/2),direction_vector[0]*np.sin(theta_rotate/2),direction_vector[1]*np.sin(theta_rotate/2),direction_vector[2]*np.sin(theta_rotate/2)]
    #print('direction_vector',direction_vector)
    #print('theta_rotate',theta_rotate)
    #print('quaternion_para',quaternion_para)
    homo = np.linalg.inv(quater2homo(quaternion_para))
    #print('homo',homo)
    return np.dot(homo,start_point - p1) + p1

                
    #after derive the H1 that from base coordinate to laser coordinate,we need to calculate the quaternion and control robot to move   
rotate_interval = 3
for i in range(120):
    center1 = np.array((253.11803385,-498.0047026,9.49255484))
    center2 = center1 + (0,0,500)
    start_point = center2 + (50,50,0)    
    rotate_point = CB_path_planning(center1,center2,start_point,np.deg2rad(i*rotate_interval))
    vector = rotate_point - center2
