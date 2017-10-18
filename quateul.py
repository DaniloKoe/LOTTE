# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 09:24:30 2017
Euler Quaternion DCM Conversion
@author: Danilo Költzsch
"""


import numpy as np


class eulerVSquat:

    
    
    def __init__(self):
        self._quat  = np.array([ 0, 0, 0, 0])
        self._euler = np.array([ 0, 0, 0])
        self._dcm   = np.array([[1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])
        # calc of euler by dcm
        self._euler = self.dirDcmToEuler(self._dcm) 
        # cal of quaternion by dcm
        self._quat = self.dirDcmToQuat(self._dcm)
        # no further commands
        
    # specifies DCM     
    def setDCM (self, dcmloc):
        self._dcm[:,:] = dcmloc[:,:]
        
    # specifies Quaternion
    def setQuat (self, quatloc):
        self._quat[:] = quatloc[:]
    
    # specifies Euler    
    def setEuler (self, euloc):
        self._euler[:] = euloc[:]
    
    # interior calc from dcm to euler
    def DcmToEuler (self):
        self._euler = self.dirDcmToEuler(self._dcm)
        
    # exterior calc from dcm to euler  
    def dirDcmToEuler (self, dcm):
        ey = np.arctan2(-dcm[0,2],dcm[2,2])
        ex = np.arcsin(dcm[1,2])
        ez = np.arctan2(-dcm[1,0],dcm[1,1])
        euler = np.array([ex, ey, ez])    
        return euler 
        
    # interior calc from dcm to quaternion      
    def DcmToQuat (self):
        q = self.dirDcmToQuat(self._dcm)
        self._quat = q
        
     # interior calc from dcm to quaternion      
    def dirDcmToQuat (self, dcm):
        cq = np.zeros(4)
        q = np.zeros(4)
        cq[0]  = (1 + dcm[0,0] + dcm[1,1] + dcm[2,2]) / 4
        cq[1]  = (1 + dcm[0,0] - dcm[1,1] - dcm[2,2]) / 4
        cq[2]  = (1 - dcm[0,0] + dcm[1,1] - dcm[2,2]) / 4
        cq[3]  = (1 - dcm[0,0] - dcm[1,1] + dcm[2,2]) / 4
         # chooses maximim element
        check = np.argmax(cq)
        if check == 0:
            q[0] = np.sqrt(cq[0])
            q[1] = ( dcm[2,1] - dcm[1,2] ) / ( 4 * q[0] )
            q[2] = ( dcm[0,2] - dcm[2,0] ) / ( 4 * q[0] )
            q[3] = ( dcm[1,0] - dcm[0,1] ) / ( 4 * q[0] )
        elif check == 1:
            q[1] = np.sqrt(cq[1])
            q[0] = ( dcm[2,1] - dcm[1,2] ) / ( 4 * q[1] )
            q[2] = ( dcm[1,0] + dcm[0,1] ) / ( 4 * q[1] )
            q[3] = ( dcm[0,2] + dcm[2,0] ) / ( 4 * q[1] )
        elif check == 2:
            q[2] = np.sqrt(cq[2])
            q[0] = ( dcm[0,2] - dcm[2,0] ) / ( 4 * q[2] )
            q[1] = ( dcm[1,0] + dcm[0,1] ) / ( 4 * q[2] )
            q[3] = ( dcm[2,1] + dcm[1,2] ) / ( 4 * q[2] )
        elif check == 3:
            q[3] = np.sqrt(cq[3])
            q[0] = ( dcm[1,0] - dcm[0,1] ) / ( 4 * q[3] )
            q[1] = ( dcm[0,2] + dcm[2,0] ) / ( 4 * q[3] )
            q[2] = ( dcm[2,1] + dcm[1,2] ) / ( 4 * q[3] )
        else:
            print('Error in Q Calc')
        return q
        
        
        
    # interior calc from quaternion to dcm   
    def QuatToDcm (self):
        self._dcm = self.dirQuatToDcm(self._quat)
    
    
    # exterior calc from quaternion to dcm   
    def dirQuatToDcm (self, quat):
        t11 = quat[0]**2 + quat[1]**2 - quat[2]**2 - quat[3]**2
        t22 = quat[0]**2 - quat[1]**2 + quat[2]**2 - quat[3]**2
        t33 = quat[0]**2 - quat[1]**2 - quat[2]**2 + quat[3]**2
        t12 = 2 * (quat[1] *  quat[2] - quat[3] * quat[0])
        t13 = 2 * (quat[1] *  quat[3] + quat[2] * quat[0])
        t21 = 2 * (quat[1] *  quat[2] + quat[3] * quat[0])
        t23 = 2 * (quat[2] *  quat[3] - quat[1] * quat[0])
        t31 = 2 * (quat[1] *  quat[3] - quat[2] * quat[0])
        t32 = 2 * (quat[2] *  quat[3] + quat[1] * quat[0])
        dcm = np.array([[t11, t12, t13],
                        [t21, t22, t23],
                        [t31, t32, t33]])
        return dcm
    

    # interior calc from euler to dcm   
    def EulerToDcm (self):
        self._dcm = self.dirEulerToDcm(self._euler)

    
    # exterior calc from euler to dcm   defined after z,y',x'' (mobil technics convention)
    def dirEulerToDcm (self, euler):
        t11 =   np.cos(euler[1]) * np.cos(euler[2]) 
        t12 =   np.cos(euler[1]) * np.sin(euler[2]) 
        t13 = - np.sin(euler[1])
        t21 =   np.sin(euler[0]) * np.sin(euler[1]) * np.cos(euler[2]) - np.cos(euler[0]) * np.sin(euler[2]) 
        t22 =   np.sin(euler[0]) * np.sin(euler[1]) * np.sin(euler[2]) + np.cos(euler[0]) * np.cos(euler[2])
        t23 =   np.sin(euler[0]) * np.cos(euler[1]) 
        t31 =   np.cos(euler[0]) * np.sin(euler[1]) * np.cos(euler[2]) + np.sin(euler[0]) * np.sin(euler[2]) 
        t32 =   np.cos(euler[0]) * np.sin(euler[1]) * np.sin(euler[2]) - np.sin(euler[0]) * np.cos(euler[2]) 
        t33 =   np.cos(euler[0]) * np.cos(euler[1])
        dcm = np.array([[t11, t12, t13],
                        [t21, t22, t23],
                        [t31, t32, t33]])
        return dcm
    
    # interior calc from quaternion to euler using previous but shortent
    def EulerToQuat (self):
        self._quat = self.dirEulerToQuat(self._euler)
        
    # exterior calc from quaternion to euler using previous but shortent
    def dirEulerToQuat (self, euler):
        dcm = self.dirEulerToDcm(euler)
        quat =  self.dirDcmToQuat(dcm)
        return quat
        
        
    # interior calc from euler to quaternion using previous but shortent 
    def QuatToEuler (self):
        self._euler = self.dirQuattoEuler(self._quat)
        
        
    # exterior calc from euler to quaternion using previous but shortent  for performence
    def dirQuatToEuler (self, quat):
        t22 = quat[0]**2 - quat[1]**2 + quat[2]**2 - quat[3]**2
        t33 = quat[0]**2 - quat[1]**2 - quat[2]**2 + quat[3]**2
        t13 = 2 * (quat[1] *  quat[3] + quat[2] * quat[0])
        t21 = 2 * (quat[1] *  quat[2] + quat[3] * quat[0])
        t23 = 2 * (quat[2] *  quat[3] - quat[1] * quat[0])
        ey = np.arctan2(-t13,t33)
        ex = np.arcsin(t23)
        ez = np.arctan2(-t21,t22)
        euler = np.array([ex, ey, ez])
        return euler
     
        
    
    # Rotation of a System defined By Direction Cosinus Matrix and By Local Euler Angles
    def SatRotFree (self, euloc, dcmloc):
        # transforms euler into quaternion
        qr2 =  self.dirEulerToQuat(euloc)
        qr  = qr2 / np.linalg.norm(qr2)
        qa2 = self.dirDcmToQuat(dcmloc)
        qa  = qa2 / np.linalg.norm(qa2)
        qak = np.array([qa[0], -qa[1], -qa[2], -qa[3]])
        # rotation of rotation into global system
        qra =   self.quatmult(qa, self.quatmult(qr, qak))
        # interior rotation for shortent text
        dcmb = self.QuatRotOfDcm(qra, dcmloc)
        angle, angleloc = self.getAngleByDcm(dcmloc, dcmb)
        # specifies global rotation rate
        #angle = self.dirQuatToEuler(qra)
        return dcmb, angle, angleloc
    
    # Rotation of a System defined By Direction Cosinus Matrix (dcmloc) By Euler Angle (rotr), Fix about a defined Axis by dcm(z Component) 
    def SatRotFix (self, rotangle, dcmloc, rotationM, vfix):
        # shortent notation by interior class
        dcmctemp = self.RotAboutFix(rotangle, dcmloc, rotationM)
        # now the magic trick by adding both worlds into one
        # chooses the fixed axis and orthogonalise the matrix to it
        case = np.array([rotationM[0,2],rotationM[1,2],rotationM[2,2]])
        check  = np.argmax(np.absolute(case))
        orient = case[check] / np.absolute(case[check])
        vfixor = vfix * orient
        vtemp = np.array([dcmctemp[0, check],dcmctemp[1, check],dcmctemp[2, check]])
        # calculates rotation by this two vectors
        qr   = self.getQuatByVect(vtemp, vfixor)
        # interior rotation for shortent text
        dcmc = self.QuatRotOfDcm(qr, dcmctemp)
        # calculating overall rotationrate
        angle, angleloc = self.getAngleByDcm(dcmloc, dcmc)
        return dcmc, angle, angleloc
    
    # Fixed to a Vector Pitch about a defined Axis (z Component)
    def SatFixed (self, dcmloc, rotationM, vfix):
        # chooses the fixed axis and orthogonalise the matrix to it
        case = np.array([rotationM[0,2],rotationM[1,2],rotationM[2,2]])
        check  = np.argmax(np.absolute(case))
        orient = case[check] / np.absolute(case[check])
        vfixor = vfix * orient
        vtemp = np.array([dcmloc[0, check],dcmloc[1, check],dcmloc[2, check]])
        # calculates rotation by this two vectors
        qr   = self.getQuatByVect(vtemp, vfixor)
        # interior rotation for shortent text
        dcmc = self.QuatRotOfDcm(qr, dcmloc)
        # calculating overall rotationrate
        angle, angleloc = self.getAngleByDcm(dcmloc, dcmc)
        return dcmc, angle, angleloc
    
    
    # Rotation about the fixed Yaw Axis for Simplification
    def RotAboutFix (self, rotangle, dcmloc, rotationM):
        # transforms euler rotation rate into quaternion
        euloc = np.array([0., 0., rotangle])
        q12 =  self.dirEulerToQuat(euloc)
        q1  =  q12 / np.linalg.norm(q12)
        # creates rotation rate about fixed axis by quaternion rotation
        qrotM2 = self.dirDcmToQuat(rotationM)
        qrotM  = qrotM2 / np.linalg.norm(qrotM2)
        qrotMk = np.array([qrotM[0], -qrotM[1], -qrotM[2], -qrotM[3]])
        #print(qrotMk)
        qr2 = self.quatmult(qrotM, self.quatmult(q1, qrotMk))
        #print(qr2)
        # doing the standard transformation like for free rot
        # transforms dcm rotation rate into quaternion
        qr  = qr2 / np.linalg.norm(qr2)
        qa2 = self.dirDcmToQuat(dcmloc)
        qa  = qa2 / np.linalg.norm(qa2)
        qak = np.array([qa[0], -qa[1], -qa[2], -qa[3]])
        # rotation of rotation into global system
        qra =   self.quatmult(qa, self.quatmult(qr, qak)) 
        # interior rotation for shortent text
        dcmctemp = self.QuatRotOfDcm(qra, dcmloc)
        return dcmctemp
    
        # Rotation about Pitch defined by Local Matrix, Transformation and Pitch Angle
    def RotAboutPitch (self, rotangle, dcmloc, rotationM):
        # transforms euler rotation rate into quaternion
        euloc = np.array([0., rotangle, 0.])
        q12 =  self.dirEulerToQuat(euloc)
        q1  =  q12 / np.linalg.norm(q12)
        # creates rotation rate about fixed axis by quaternion rotation
        qrotM2 = self.dirDcmToQuat(rotationM)
        qrotM  = qrotM2 / np.linalg.norm(qrotM2)
        qrotMk = np.array([qrotM[0], -qrotM[1], -qrotM[2], -qrotM[3]])
        qr2 = self.quatmult(qrotM, self.quatmult(q1, qrotMk))
        # doing the standard transformation like for free rot
        # transforms dcm rotation rate into quaternion
        qr  = qr2 / np.linalg.norm(qr2)
        qa2 = self.dirDcmToQuat(dcmloc)
        qa  = qa2 / np.linalg.norm(qa2)
        qak = np.array([qa[0], -qa[1], -qa[2], -qa[3]])
        # rotation of rotation into global system
        qra =   self.quatmult(qa, self.quatmult(qr, qak))
        qrak = np.array([qra[0], -qra[1], -qra[2], -qra[3]]) 
        # arrange vector as quaternion
        v1 = np.array([0, dcmloc[0,0], dcmloc[1,0], dcmloc[2,0]])
        v2 = np.array([0, dcmloc[0,1], dcmloc[1,1], dcmloc[2,1]])
        v3 = np.array([0, dcmloc[0,2], dcmloc[1,2], dcmloc[2,2]])
        # rotation of the dcm vector as quaternion
        q1   = self.quatmult(qra, self.quatmult(v1, qrak))
        q2   = self.quatmult(qra, self.quatmult(v2, qrak))
        q3   = self.quatmult(qra, self.quatmult(v3, qrak))  
        q1n  = q1 / np.linalg.norm(q1)
        q2n  = q2 / np.linalg.norm(q2)
        q3n  = q3 / np.linalg.norm(q3)
        # transform quaternion into dcm
        dcmctemp = np.array([[q1n[1], q2n[1], q3n[1]],
                             [q1n[2], q2n[2], q3n[2]],
                             [q1n[3], q2n[3], q3n[3]]])
        return dcmctemp
    
    
        # Rotation of DCM as Vector Components defined by global Quaternion
    def QuatRotOfDcm(self, qr, dcma):
        qrn  = qr / np.linalg.norm(qr)
        qrnk = np.array([qrn[0], -qrn[1], -qrn[2], -qrn[3]])
        # arrange vector as quaternion
        v1 = np.array([0, dcma[0,0], dcma[1,0], dcma[2,0]])
        v2 = np.array([0, dcma[0,1], dcma[1,1], dcma[2,1]])
        v3 = np.array([0, dcma[0,2], dcma[1,2], dcma[2,2]])
        # rotation of the dcm vector as quaternion
        q1   = self.quatmult(qrn, self.quatmult(v1, qrnk))
        q2   = self.quatmult(qrn, self.quatmult(v2, qrnk))
        q3   = self.quatmult(qrn, self.quatmult(v3, qrnk))  
        q1n  = q1 / np.linalg.norm(q1)
        q2n  = q2 / np.linalg.norm(q2)
        q3n  = q3 / np.linalg.norm(q3)
        # transform quaternion into dcm
        dcmctemp = np.array([[q1n[1], q2n[1], q3n[1]],
                             [q1n[2], q2n[2], q3n[2]],
                             [q1n[3], q2n[3], q3n[3]]])
        return dcmctemp
    
    # Rotation of DCM as Vector Components defined by global Quaternion
    def QuatRotOfV(self, qr, va):
        qrn  = qr / np.linalg.norm(qr)
        qrnk = np.array([qrn[0], -qrn[1], -qrn[2], -qrn[3]])
        # arrange vector as quaternion
        v1 = np.array([0, va[0], va[1], va[2]])
        # rotation of the dcm vector as quaternion
        q1   = self.quatmult(qrn, self.quatmult(v1, qrnk))
        q1n  = q1 / np.linalg.norm(q1)
        # transform quaternion into dcm
        vb   = np.array([q1n[1], q1n[2], q1n[3]])
        return vb

    
        # Rotates DCM Vector defined By DCM
    def DcmRotOfDcm(self, dcmrot, dcma):
        qr   = self.dirDcmToQuat(dcmrot)
        dcmb = self.QuatRotOfDcm(qr, dcma)
        return dcmb
    
    
        # Rotates Euler Rotation defined By DCM
    def DcmRotOfEuler(self, dcmrot, eulerA):
        qa     = self.dirEulerToQuat(eulerA)
        qr2    = self.dirDcmToQuat(dcmrot)
        qr     = qr2 / np.linalg.norm(qr2)
        qrk    = np.array([qr[0], -qr[1], -qr[2], -qr[3]])
        qb     = self.quatmult(qr, self.quatmult(qa, qrk))
        eulerB = self.dirQuatToEuler(qb)
        return eulerB
         # Rotates DCM Vector defined By DCM
    
        # Reverse Rotation of Vector defined By DCM     
    def DcmRotOfV(self, dcmrot, va):
        qr   = self.dirDcmToQuat(dcmrot)
        vb = self.QuatRotOfV(qr, va)
        return vb  
    
         # Reverse Rotation of  DCM Vector defined By DCM
    def revDcmRotOfDcm(self, dcmrot, dcma):
        qr2  = self.dirDcmToQuat(dcmrot)
        qr   = qr2 / np.linalg.norm(qr2)
        qrk   = np.array([qr[0], -qr[1], -qr[2], -qr[3]])
        dcmb = self.QuatRotOfDcm(qrk, dcma)
        return dcmb
    
    
        # Reverse Rotation of Euler Rotation defined By DCM
    def revDcmRotOfEuler(self, dcmrot, eulerA):
        qa     = self.dirEulerToQuat(eulerA)
        qr2    = self.dirDcmToQuat(dcmrot)
        qr     = qr2 / np.linalg.norm(qr2)
        qrk    = np.array([qr[0], -qr[1], -qr[2], -qr[3]])
        qb     = self.quatmult(qrk, self.quatmult(qa, qr))
        eulerB = self.dirQuatToEuler(qb)
        return eulerB
    
     # Reverse Rotation of Vector defined By DCM     
    def revDcmRotOfV(self, dcmrot, va):
        qr   = self.dirDcmToQuat(dcmrot)
        qrk  = np.array([qr[0], -qr[1], -qr[2], -qr[3]])
        vb = self.QuatRotOfDcm(qrk, va)
        return vb  
    
    def getQuatByVect(self, v1, v2):
        v3 = v2 - v1
        # building dcm matrix by two vectors and transform into q
        v1dcm = self.buildDcmByV(v1, v3)
        q12 = self.dirDcmToQuat(v1dcm)
        q1  = q12 / np.linalg.norm(q12)
        q1k = np.array([q1[0], -q1[1], -q1[2], -q1[3]]) 
        vq2 = np.array([0, v2[0], v2[1], v2[2]])
        # rotation q into local syste, which specifies rotation rate about y
        q2  = self.quatmult(q1k, self.quatmult(vq2, q1))
        # shall be pitch angle 0
        vz2 = np.array([q2[1], 0, q2[3]])
        vz  = vz2 / np.linalg.norm(vz2)
        vy  = np.array([0, 1, 0])
        # rebuild x
        vx  = np.cross(vy,vz)
        # rebuild dcmr in local system
        dcmr = np.array([[vx[0],vy[0],vz[0]],
                         [vx[1],vy[1],vz[1]],
                         [vx[2],vy[2],vz[2]]])
        qrloc = self.dirDcmToQuat(dcmr)
        # back to global
        qr = self.quatmult(q1, self.quatmult(qrloc, q1k))
        return qr
        
    
    # local Rotoation of Quaternion    
    def quatmult (self, q1, q2):
        qs = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3]
        qx = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2]
        qy = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1]
        qz = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0]
        q = np.array([qs, qx, qy, qz])
        return q
    
    # Rotation between too systems by dcm in  global and local coordiantes
    def getAngleByDcm (self, dcma, dcmb):
        # calculation of rotation forced by fixed axis
        dcma2 = np.matrix([[dcma[0,0], dcma[0,1], dcma[0,2]],
                           [dcma[1,0], dcma[1,1], dcma[1,2]],
                           [dcma[2,0], dcma[2,1], dcma[2,2]]])
        dcmb2 = np.matrix([[dcmb[0,0], dcmb[0,1], dcmb[0,2]],
                           [dcmb[1,0], dcmb[1,1], dcmb[1,2]],
                           [dcmb[2,0], dcmb[2,1], dcmb[2,2]]])
        dcmfixrot = dcmb2 * dcma2.transpose()
        # forced fixed rotation into quaternion in global coords
        qrfix2  = self.dirDcmToQuat(dcmfixrot)
        qrfix   = qrfix2 / np.linalg.norm(qrfix2)
        # actual position of satellite into local system
        qloc2 = self.dirDcmToQuat(dcma)  
        qloc  = qloc2 / np.linalg.norm(qloc2)
        qlock = np.array([qloc[0], -qloc[1], -qloc[2], -qloc[3]])
        # specifies local and global rotation rate
        qrloc    = self.quatmult(qlock, self.quatmult(qrfix, qloc)) 
        angle    = self.dirQuatToEuler(qrfix)
        angleloc = self.dirQuatToEuler(qrloc)
        return angle, angleloc
    
        # Return Angle between Vectors in Deg
    def getAngleByV(self, va, vb):
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        if cos_angle >= 1:
            angle = 0.
        elif cos_angle <= -1:
            angle = 180.
        else:
            angle = np.arccos(cos_angle) * 180 / np.pi

        return angle
    
            # Return Angle between Vectors in Rad
    def getAngleRadByV(self, va, vb):
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        if cos_angle >= 1:
            angle = 0.
        elif cos_angle <= -1:
            angle = np.pi / 2
        else:
            angle       = np.arccos(cos_angle)
        return angle
    
         # Return Cosinus of Angle between Vectors
    def getCosByV(self, va, vb):
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        return cos_angle 
    
    
         # Return Cosinus of Sun Angle between Vectors
    def getRadiCosByV(self, va, vb):
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        if cos_angle < 0:
            cos_angle = 0
        return cos_angle 
    
            # Return Angle between Vectors in Deg (second Vector by z Component of DCM)
    def getAngleByVZofDcm(self, va, dcmb):
        vb          = np.array([dcmb[0,2],dcmb[1,2],dcmb[2,2]])
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        angle       = np.arccos(cos_angle) * 180 / np.pi
        return angle
    
         # Return Cosinus of Angle between Vectors (second Vector by z Component of DCM)
    def getCosByVZofDcm(self, va, dcmb):
        vb          = np.array([dcmb[0,2],dcmb[1,2],dcmb[2,2]])
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        return cos_angle 
    
        
         # Return Cosinus of incomming Radiation (upto 90°), Angle between Vectors (second Vector by z Component of DCM)
    def getRadiCosByVZofDcm(self, va, dcmb):
        vb          = np.array([dcmb[0,2],dcmb[1,2],dcmb[2,2]])
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        if cos_angle < 0:
            cos_angle = 0
        return cos_angle 
    
             # Return Cosinus of incomming Radiation (upto 90°), Angle for shader in two axis devided between Vectors (second Vector by z Component of DCM)
    def getRadiShadCosByVZofDcm(self, va, dcmb):
        vb          = np.array([dcmb[0,2],dcmb[1,2],dcmb[2,2]])
        dot         = np.dot(va,vb)
        va_modulus  = np.sqrt((va**2).sum())
        vb_modulus  = np.sqrt((vb**2).sum())
        cos_angle   = dot / va_modulus / vb_modulus
        if cos_angle < 0:
            cos_angle1 = 0
            cos_angle2 = 0
        elif cos_angle == 1:
            cos_angle1 = 1
            cos_angle2 = 1
        else:
            vapr1, vapr2    = self.getProjOfDcmAndV(va, dcmb)
            dot1         = np.dot(vapr1,vb)
            va_modulus1  = np.sqrt((vapr1**2).sum())
            vb_modulus1  = np.sqrt((vb**2).sum())
            cos_angle1   = dot1 / va_modulus1 / vb_modulus1
            dot2         = np.dot(vapr2,vb)
            va_modulus2  = np.sqrt((vapr2**2).sum())
            vb_modulus2  = np.sqrt((vb**2).sum())
            cos_angle2   = dot2 / va_modulus2 / vb_modulus2
        return cos_angle1, cos_angle2  
    
    
                 # Return Angle of incomming Radiation (upto 90°), Angle for shader in two axis devided between Vectors (second Vector by z Component of DCM)
    def getShaderAngleByVZofDcm(self, va, dcmb):
        vb1          = np.array([dcmb[0,0],dcmb[1,0],dcmb[2,0]])
        vb2          = np.array([dcmb[0,1],dcmb[1,1],dcmb[2,1]])
        vb3          = np.array([dcmb[0,2],dcmb[1,2],dcmb[2,2]])
        vapr1, vapr2    = self.getProjOfDcmAndV(va, dcmb)
        if sum(vapr1**2) == 0:
           angle1 = 0
        else:
           dot1         = np.dot(vapr1,vb3)
           dot12        = np.dot(vapr1,vb1)
           va_modulus1  = np.sqrt((vapr1**2).sum())
           vb_modulus1  = np.sqrt((vb3**2).sum())
           if dot12 < 0:
              angle1   = - np.arccos(dot1 / va_modulus1 / vb_modulus1)
           else:
              angle1   =   np.arccos(dot1 / va_modulus1 / vb_modulus1)
        if sum(vapr2**2) == 0:
              angle2 = 0
        else:
              dot2         = np.dot(vapr2,vb3)
              dot22        = np.dot(vapr2,vb2)
              va_modulus2  = np.sqrt((vapr2**2).sum())
              vb_modulus2  = np.sqrt((vb3**2).sum())
              if dot22 < 0:
                 angle2   = - np.arccos(dot2 / va_modulus2 / vb_modulus2)
              else:
                 angle2   =   np.arccos(dot2 / va_modulus2 / vb_modulus2)
        return angle1, angle2 
    
        # returns a projection of a vektor into a plane defined by a dcm plane1 (v1,v3) and plane2 (v2,v3)
    def getProjOfDcmAndV(self,va,dcmb):
        vb1 = dcmb[0,:]
        vb2 = dcmb[1,:]
        vb3 = dcmb[2,:]
        v1 = sum(va * vb1)
        v2 = sum(va * vb2)
        v3 = sum(va * vb3)
        u1 = sum(vb1 * vb1)
        u2 = sum(vb2 * vb2)
        u3 = sum(vb3 * vb3)
        vpr1 = v1 / u1 * vb1 + v3 / u3 * vb3
        vpr2 = v2 / u2 * vb2 + v3 / u3 * vb3
        return vpr1, vpr2
    
    
    
    # Rotation between too systems by dcm in  global and local coordiantes in quaternion
    def getQuatByDcm (self, dcma, dcmb):
        # calculation of rotation forced by fixed axis
        dcma2 = np.matrix([[dcma[0,0], dcma[0,1], dcma[0,2]],
                           [dcma[1,0], dcma[1,1], dcma[1,2]],
                           [dcma[2,0], dcma[2,1], dcma[2,2]]])
        dcmb2 = np.matrix([[dcmb[0,0], dcmb[0,1], dcmb[0,2]],
                           [dcmb[1,0], dcmb[1,1], dcmb[1,2]],
                           [dcmb[2,0], dcmb[2,1], dcmb[2,2]]])
        dcmfixrot = dcmb2 * dcma2.transpose()
        # forced fixed rotation into quaternion in global coords
        qrfix2  = self.dirDcmToQuat(dcmfixrot)
        qrfix   = qrfix2 / np.linalg.norm(qrfix2)
        # actual position of satellite into local system
        qloc2 = self.dirDcmToQuat(dcma)  
        qloc  = qloc2 / np.linalg.norm(qloc2)
        qlock = np.array([qloc[0], -qloc[1], -qloc[2], -qloc[3]])
        # specifies local and global rotation rate
        qrloc    = self.quatmult(qlock, self.quatmult(qrfix, qloc)) 
        return qrfix, qrloc
    
    
    # Rotation between too systems by quaternion in global and local coordiantes
    def getAngleByQuat(self, quata, quatb):
        # transform quat to dcm
        qa = quata / np.linalg.norm(quata)
        qb = quatb / np.linalg.norm(quatb)
        dcma = self.dirQuatToDcm(qa)
        dcmb = self.dirQuatToDcm(qb)
        # calculation of rotation forced by fixed axis
        dcma2 = np.matrix([[dcma[0,0], dcma[0,1], dcma[0,2]],
                           [dcma[1,0], dcma[1,1], dcma[1,2]],
                           [dcma[2,0], dcma[2,1], dcma[2,2]]])
        dcmb2 = np.matrix([[dcmb[0,0], dcmb[0,1], dcmb[0,2]],
                           [dcmb[1,0], dcmb[1,1], dcmb[1,2]],
                           [dcmb[2,0], dcmb[2,1], dcmb[2,2]]])
        dcmrot = dcmb2 * dcma2.transpose()
        # forced fixed rotation into quaternion in global coords
        qr2  = self.dirDcmToQuat(dcmrot)
        qr = qr2 / np.linalg.norm(qr2)
        # actual position of satellite into local system
        qa = quata / np.linalg.norm(quata)
        qak = np.array([qa[0], -qa[1], -qa[2], -qa[3]])
        # specifies local and global rotation rate
        qrloc    = self.quatmult(qak, self.quatmult(qr, qa)) 
        angle    = self.dirDcmToEuler(dcmrot)
        angleloc = self.dirQuatToEuler(qrloc)
        return angle, angleloc
    
    # builds DCM by fixed and flight vector
    def buildDcmByV (self, v3, v1):
        v3n = v3 / np.linalg.norm(v3)
        v1n = v1 / np.linalg.norm(v1)
        v2n  = np.cross(v3n,v1n)
        v11n = np.cross(v2n,v3n)
                # has to be multiplied on the right site
        dcm = np.array([[v11n[0],   v2n[0],  v3n[0]],
                        [v11n[1],   v2n[1],  v3n[1]],
                        [v11n[2],   v2n[2],  v3n[2]]])
        return dcm
    
    # norms DCM
    def correctDcm (self, dcma):
        # arrange dcm as vector
        v1 = np.array([dcma[0,0], dcma[1,0], dcma[2,0]])
        v3 = np.array([dcma[0,2], dcma[1,2], dcma[2,2]])
        # correct dcm by interior class
        dcmb = self.buildDcmByV(v3, v1)
        return dcmb
    
  # builds Normal vector by thre points
    def buildNormalByPoint(self, p3, p2, p1):
        v1  = p2 - p1
        v3  = p2 - p3
        v3n = v3 / np.linalg.norm(v3)
        v1n = v1 / np.linalg.norm(v1)
        vn  = np.cross(v3n,v1n)
        return vn
        