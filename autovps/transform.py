import numpy as np


class Transform(object):
    """Transform matrix class
    """
    tform = np.eye(4)

    def __init__(self, tform):
        self.tform = np.asarray(tform)

    def get_matrix(self):
        return(self.tform)

    def get_inverse(self):
        return Transform(np.linalg.pinv(self.get_matrix()))
        
    def rotx(self, theta):
        """ Apply a rotation of theta degrees around the x axis
        """

        theta = np.radians(theta)
        R = np.asarray([[1, 0, 0, 0],
                        [0, np.cos(theta), -np.sin(theta), 0],
                        [0, np.sin(theta), np.cos(theta), 0],
                        [0, 0, 0, 1]])

        self.concatenate_matrix(R)

    def roty(self, theta):
        """ Apply a rotation of theta degrees around the y axis
        """
        theta = np.radians(theta)
        R = np.asarray([[np.cos(theta), 0, np.sin(theta), 0],
                        [0, 1, 0, 0],
                        [-np.sin(theta), 0, np.cos(theta), 0],
                        [0, 0, 0, 1]] )

        self.concatenate_matrix(R)

    def rotz(self, theta):
        """ Apply a rotation of theta degrees around the z axis
        """
        theta = np.radians(theta)
        R = np.asarray([[np.cos(theta), -np.sin(theta), 0, 0],
                        [np.sin(theta), np.cos(theta), 0, 0],
                        [0, 0, 1, 0],
                        [0, 0, 0, 1]])

        self.concatenate_matrix(R)

    def scaleXYZ(self, x, y, z):
        """Apply a scaling given by x,y,z
        """
        S = np.diag([x, y, z, 1])
        self.concatenate_matrix(S)


    def scale(self, s):
        """Apply a scaling given by the 1x3 vector s
        """
        self.scaleXYZ(s[0], s[1], s[2])


    def get_scale(self):
        dims = [np.sqrt(np.dot(self.tform[:, i].T.tolist(), self.tform[:, i].tolist())) for i in range(3)]
        return(np.squeeze(dims))

    def get_position(self):
        return(np.asarray(self.tform[0:3, 3]).squeeze())

    def concatenate_matrix(self, tform):
        """Right multiply the current transform by tform

        Args:
            tform (:obj:`numpy.ndarray`): a 4x4 transform to apply

        """

        self.tform = np.dot(self.tform, tform)


    def translateXYZ(self, x, y, z):
        self.tform[0:3, 3] = self.tform[0:3, 3]+np.matrix([[x, y, z]]).T

    # overload operators, no copy
    def __mul__(self, tform):
        T2 = Transform(self.get_matrix())
        T2.concatenate_matrix(tform.get_matrix())
        return(T2)

    def get_rotation(self):
        R = Transform(self.get_matrix())
        R.scale(1.0/self.get_scale())
        return(R.tform[0:3, 0:3])

    def siemens_orientation(self):

        """Given a 3-item list (or other iterable) that represents a normal vector
        to the "imaging" plane, this function determines the orientation of the
        vector in 3-dimensional space. It returns a tuple of (angle, orientation)
        in which angle is e.g. "Tra" or "Tra>Cor -6" or "Tra>Sag 14.1 >Cor 9.3"
        and orientation is e.g. "Sag" or "Cor-Tra".

        For double angulation, errors in secondary angle occur that may be due to
        rounding errors in internal Siemens software, which calculates row and
        column vectors.
        """
        # from hcpre

        # docstring paraphrases IDL comments
        TOLERANCE = 1.e-4
        orientations = ('Sag', 'Cor', 'Tra')

        final_angle = ""

        # [IDL] evaluate orientation of normal vector:
        #
        # Find principal direction of normal vector (i.e. axis with its largest
        # component)
        # Find secondary direction (second largest component)
        # Calc. angle btw. projection of normal vector into the plane that
        #     includes both principal and secondary directions on the one hand
        #     and the principal direction on the other hand ==> 1st angulation:
        #     "principal>secondary = angle"
        # Calc. angle btw. projection into plane perpendicular to principal
        #     direction on the one hand and secondary direction on the other
        #     hand ==> 2nd angulation: "secondary>third dir. = angle"

        # get principal, secondary and ternary directions

        normal = np.squeeze(np.array(self.tform[0:3, 2]))
        ternary, secondary, principal = np.argsort(np.abs(normal))

        # [IDL] calc. angle between projection into third plane (spawned by
        # principle & secondary directions) and principal direction:
        angle_1 = np.degrees(np.arctan2(normal[secondary], normal[principal]))

        # [IDL] calc. angle btw. projection on rotated principle direction and
        # secondary direction:
        # projection on rotated principle dir.
        new_normal_ip = np.sqrt((normal[principal] ** 2) + (normal[secondary] ** 2))

        angle_2 = np.degrees(np.arctan2(normal[ternary], new_normal_ip))

        # [IDL] SIEMENS notation requires modifications IF principal dir. indxs SAG !
        # [PS] In IDL, indxs is the name of the variable that is "secondary" here.
        #      Even with that substitution, I don't understand the comment above.
        if not principal:
            if abs(angle_1) > 0:
                sign1 = angle_1 / abs(angle_1)
            else:
                sign1 = 1.0

            angle_1 -= (sign1 * 180.0)
            angle_2 *= -1

        if (abs(angle_2) < TOLERANCE) or (abs(abs(angle_2) - 180) < TOLERANCE):
            if (abs(angle_1) < TOLERANCE) or (abs(abs(angle_1) - 180) < TOLERANCE):
                # [IDL] NON-OBLIQUE:
                final_angle = orientations[principal]
                final_orientation = ""
            else:
                # [IDL] SINGLE-OBLIQUE:
                final_angle = "%s>%s %.1f" % \
                        (orientations[principal], orientations[secondary],
                         (1 * angle_1)
                        )
                final_orientation = orientations[principal] + '-' + orientations[secondary]
        else:
            # [IDL] DOUBLE-OBLIQUE:
            final_angle = "%s>%s %.1f >%s %.1f" % \
                    (orientations[principal], orientations[secondary],
                     (1 * angle_1), orientations[ternary], (1 * angle_2))
            final_orientation = "%s-%s-%s" % \
                    (orientations[principal], orientations[secondary],
                     orientations[ternary])

        # inplane rotation
        R = self.get_rotation()

        if orientations[principal] == 'Tra':
            ip_rotation = np.degrees(np.arctan2(-R[0, 1], R[0, 0]))
        if orientations[principal] == 'Cor':
            ip_rotation = np.degrees(np.arctan2(R[2, 1], -R[2, 0]))
        if orientations[principal] == 'Sag':
            ip_rotation = 90.0 - np.degrees(np.arctan2(-R[1, 1], -R[1, 0]))

        return((final_angle, ip_rotation))
