# Monte-Carlo photon transport code for cylindrical samples
import scipy.sparse as sp
import sys, os
import numpy as np
import random as rand
import math

class MC_light_absorption:
    # Constructor
    def __init__(self,InputLines):
        # Create input dictionary
        input_dict = {}
        
        for l in InputLines:
            thisSet = l.split("-")
            thisSet[0] = thisSet[0].rstrip(' ')
            thisSet[1] = thisSet[1].rstrip('\n').lstrip(' ').rstrip('\r')
            
            input_dict[thisSet[0]] = thisSet[1]
            
        # Now get the values from dictionary
        self.mu_s = float(input_dict['Scattering coefficient'])
        self.mu_c = float(input_dict['Chemical absorption coefficient'])
        self.mu_b = float(input_dict['Background absorption coefficient'])
        self.g = float(input_dict['Scattering anisotropy'])
        self.eta = float(input_dict['Refractive index'])
        self.photonSet = int(input_dict['No. of photons'])
        self.tSteps = int(input_dict['Total number of steps'])
        self.wTh = float(input_dict['Weight threshold'])
        self.m = int(input_dict['Roulette constant'])
        self.output_name = input_dict['Output casename']
        
        self.sample_h = float(input_dict['Sample height'])
        self.sample_r = float(input_dict['Sample radius'])
        self.eta_0 = float(input_dict['Refractive index (at Z=0)'])
        self.eta_h = float(input_dict['Refractive index (at Z=height)'])
        self.eta_r = float(input_dict['Refractive index (at X^2 + Y^2 = R^2)'])
        self.n_segs = int(input_dict['No. of segments'])
        self.dh = self.sample_h/self.n_segs
        
    def initializeStructs(self):
        # Total attenuation coefficient
        self.mu_t = self.mu_b + self.mu_c + self.mu_s
        self.mu_a = self.mu_b + self.mu_c
        # Total photon source weight
        self.total_source_wt = 0.0
        # Total number of photons
        self.N_pt = self.photonSet
       
        # Absorption field
        self.absorbed_field = np.zeros(self.n_segs)
        
        self.beer_lambert_absorption = 1.0 - math.exp(-self.mu_c*self.sample_h)
        
        print self.beer_lambert_absorption
    
    def timeLoop(self):
        self.releasedPhoton = 0.0
        while self.N_pt>0:
            coords, cosines = self.select_starting_pt()
            self.N_pt += -1
            
            current_wt = 1.0
            
            reflections = 0
            
            while current_wt>0.0:
                step = self.sample_step()
                
                trial_coords, position = self.move_photon(coords,cosines,step)
                
                if position=='inside':
                    coords = list(trial_coords)
                    self.absorb_wts(trial_coords,current_wt)
                    cosines = self.update_cosines(cosines)
                    current_wt *= (1.0 - self.mu_a/self.mu_t)
                else:
                    while position!='inside' and current_wt!=0 and reflections<100:
                        coords, trial_coords, step, cosines, current_wt, position = self.reflect_photon(coords,trial_coords,step,cosines,current_wt,position)
                        reflections += 1

                    if position=='inside' and current_wt!=0.0:
                        coords = list(trial_coords)
                        self.absorb_wts(trial_coords,current_wt)
                        cosines = self.update_cosines(cosines)
                        current_wt *= (1.0 - self.mu_a/self.mu_t)
                        
                if current_wt<self.wTh:
                    if rand.uniform(0,1.0)<1.0/float(self.m):
                        current_wt = 0
                
                if current_wt==0:
                    print 'Set ', self.N_pt, ' completed!'
                    break
        
    def select_starting_pt(self):
        xi_1 = rand.uniform(0,self.sample_r)
        xi_2 = rand.uniform(0,2*math.pi)
        
        coords = np.zeros(3)
        cosines = np.zeros(3)
        
        coords[0] = xi_1*math.cos(xi_2)
        coords[1] = xi_1*math.sin(xi_2)
        cosines[2] = 1.0
        
        return coords, cosines
    
    def move_photon(self,coords,cosines,step):
        trial = np.zeros(3)
        
        trial = coords + step*cosines
        
        position = self.check_boundaries(trial)

        return trial, position            
    
    def sample_step(self):
        xi = rand.uniform(1e-16,1.0)
        
        step = -math.log(xi)/self.mu_t
        
        return step
    
    def check_boundaries(self,trial_coords):
        if trial_coords[2]<0.0:
            position = 'boundary1'
        elif trial_coords[2]>self.sample_h:
            position = 'boundary2'
        elif (trial_coords[0]**2 + trial_coords[1]**2)>self.sample_r**2:
            position = 'boundary3'
        else:
            position = 'inside'
            
        return position
    
    def update_cosines(self,cosines):
        xi_1 = rand.uniform(0.0,1.0)
        
        cosT = (1.0/(2.0*self.g))*(1+self.g**2 - ((1-self.g**2)/(1-self.g+2*self.g*xi_1))**2);
        sinT = math.sqrt(1.0 - cosT**2)
        
        phi = rand.uniform(0,1.0)*2*math.pi
        sinPhi, cosPhi = math.sin(phi), math.cos(phi)
        
        if cosines[2] != 1.0:
            cosines[0] = (sinT/np.sqrt(1-cosines[2]**2))*(cosines[0]*cosines[2]*cosPhi - cosines[1]*sinPhi) + cosines[0]*cosT
            cosines[1] = (sinT/np.sqrt(1-cosines[2]**2))*(cosines[1]*cosines[2]*cosPhi + cosines[0]*sinPhi) + cosines[1]*cosT
            cosines[2] = -np.sqrt(1-cosines[2]**2)*sinT*cosPhi + cosines[2]*cosT
        else:
            cosines[0] = sinT*cosPhi
            cosines[1] = sinT*sinPhi
            cosines[2] = cosT

        reFactor = np.linalg.norm(cosines,2)
        
        if reFactor > 1.0:
            cosines *= 1.0/reFactor
            
        return cosines
    
    def absorb_wts(self,coords,current_wt):
        idx = int(coords[2]/self.dh)
        
        self.absorbed_field[idx] += current_wt*(self.mu_c/self.mu_t)
        
    def reflect_photon(self,coords,trial_coords,step,cosines,current_wt,position):
        surface_coords = np.zeros(3)
        new_cosines = np.zeros(3)
        new_coords = np.zeros(3)
        normal = np.zeros(3)
        reflected_length = 0.0
        
        if position=='boundary1':
            theta = math.acos(-cosines[2])
            sinTheta_e = (self.eta/self.eta_0)*math.sin(theta)
            
            if sinTheta_e == 0.0:
                Ri = 0.0
            elif sinTheta_e < 1.0:
                theta_e = math.asin(sinTheta_e)
                Ri = (math.sin(theta - theta_e)/math.sin(theta + theta_e))**2
                Ri = 0.5*(Ri + (math.tan(theta - theta_e)/math.tan(theta + theta_e))**2)        
            else:
                Ri = 1.0
                
            if rand.uniform(0,1.0)>Ri:
                current_wt = 0.0
            else:
                reflected_length = trial_coords[2]/cosines[2]
                
                surface_coords[0] = trial_coords[0] - reflected_length*cosines[0]
                surface_coords[1] = trial_coords[1] - reflected_length*cosines[1]
                
                new_cosines[0], new_cosines[1] = cosines[0], cosines[1]
                new_cosines[2] = -cosines[2]
                
                new_coords[0], new_coords[1] = trial_coords[0], trial_coords[1]
                new_coords[2] = -trial_coords[2]
            
            position = self.check_boundaries(new_coords)
        elif position=='boundary2':
            theta = math.acos(cosines[2])
            sinTheta_e = (self.eta/self.eta_h)*math.sin(theta)
            
            if sinTheta_e == 0.0:
                Ri = 0.0
            elif sinTheta_e < 1.0:
                theta_e = math.asin(sinTheta_e)
                Ri = (math.sin(theta - theta_e)/math.sin(theta + theta_e))**2
                Ri = 0.5*(Ri + (math.tan(theta - theta_e)/math.tan(theta + theta_e))**2)        
            else:
                Ri = 1.0
                
            if rand.uniform(0,1.0)>Ri:
                current_wt = 0.0
            else:    
                reflected_length = (trial_coords[2]-self.sample_h)/cosines[2]
                
                surface_coords[0] = trial_coords[0] - reflected_length*cosines[0]
                surface_coords[1] = trial_coords[1] - reflected_length*cosines[1]
                
                new_cosines[0], new_cosines[1] = cosines[0], cosines[1]
                new_cosines[2] = -cosines[2]
                
                new_coords[0], new_coords[1] = trial_coords[0], trial_coords[1]
                new_coords[2] = self.sample_h - trial_coords[2]

            position = self.check_boundaries(new_coords)
        else:
            surface_coords = self.find_intersection(coords,trial_coords)
            
            reflected_length = np.linalg.norm(trial_coords-surface_coords,2)
            
            normal[0],normal[1] = surface_coords[0], surface_coords[1]
            normal *= 1.0/np.linalg.norm(normal,2)

            cosTheta = np.dot(normal,cosines)
            theta = math.acos(cosTheta)
            sinTheta_e = (self.eta/self.eta_r)*math.sin(theta)
            
            if sinTheta_e == 0.0:
                Ri = 0.0
            elif sinTheta_e < 1.0:
                theta_e = math.asin(sinTheta_e)
                Ri = (math.sin(theta - theta_e)/math.sin(theta + theta_e))**2
                Ri = 0.5*(Ri + (math.tan(theta - theta_e)/math.tan(theta + theta_e))**2)        
            else:
                Ri = 1.0
            
            if rand.uniform(0,1.0)>Ri:
                current_wt = 0.0
            else:
                new_cosines = cosines - 2.0*normal*cosTheta
                new_cosines *= 1.0/np.linalg.norm(new_cosines,2)
                new_coords = surface_coords + reflected_length*new_cosines
            
            position = self.check_boundaries(new_coords)

        return surface_coords, new_coords, reflected_length, new_cosines, current_wt, position
    
    def find_intersection(self,coords,trial_coords):
        x_1, y_1 = coords[0], coords[1]
        t_x, t_y = trial_coords[0]-coords[0], trial_coords[1]-coords[1]
        
        a = t_x**2 + t_y**2
        b = 2*(x_1*t_x + y_1*t_y)
        c = x_1**2 + y_1**2 - self.sample_r**2
        
        D = b**2 - 4*a*c
        
        alpha_1 = (-b + math.sqrt(D))/(2*a)
        alpha_2 = (-b - math.sqrt(D))/(2*a)
        
        if alpha_1>0:
            surface_pt = coords + alpha_1*(trial_coords - coords)
        else:
            surface_pt = coords - alpha_2*(trial_coords - coords)
            
        return surface_pt
    
    def write_outputs(self):
        ofile = open(self.output_name+'.csv','w')
        
        print >> ofile, 'Location,Exciation absorption'
        
        for segs in xrange(0,self.n_segs):
            print >> ofile, str(0.5*self.dh*(2*segs+1))+','+str(self.absorbed_field[segs]/(self.photonSet*self.beer_lambert_absorption))
            
        ofile.close()