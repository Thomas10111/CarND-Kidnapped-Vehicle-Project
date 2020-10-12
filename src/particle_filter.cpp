/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include "map.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
  /**
   	* TODO: Set the number of particles. Initialize all particles to 
   	*   first position (based on estimates of x, y, theta and their uncertainties
   	*   from GPS) and all weights to 1. 
   	* TODO: Add random Gaussian noise to each particle.
   	* NOTE: Consult particle_filter.h for more information about this method 
   	*   (and others in this file).
   	*/
  num_particles = 30;  // TODO: Set the number of particles
  std::default_random_engine gen;

  // Standard deviations for x, y, and theta
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  for( int i=0; i<num_particles; ++i)
  {   
    Particle p = {i, dist_x(gen), dist_y(gen), dist_theta(gen), 1, {}, {}, {}};
    particles.push_back( p );
  }
  is_initialized = true;
  
  for(Particle p: particles)
  {   
  	std::cout << "p: " << p.id << "   x: " << p.x << "   y: " << p.y << std::endl;
  }  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::default_random_engine gen;

  // Standard deviations for x, y, and theta
  const double std_x = std_pos[0];
  const double std_y = std_pos[1];
  const double std_theta = std_pos[2];
  
  for(Particle &p: particles)
  {
    double x = p.x;
    double y = p.y;
    double theta = p.theta;// rad
    double theta_dot = yaw_rate;// rad/s

    x = x + velocity/theta_dot * ( sin(theta + yaw_rate*delta_t) - sin(theta) );
    y = y + velocity/theta_dot * ( cos(theta) - cos(theta + yaw_rate*delta_t) );
    theta = theta + yaw_rate * delta_t;
    
    std::normal_distribution<double> dist_x(x, std_x);
  	std::normal_distribution<double> dist_y(y, std_y);
  	std::normal_distribution<double> dist_theta(theta, std_theta);
    
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);    
  }
   	
  for( Particle const &p: particles)
  {
    std::cout<<"p.id: " << p.id << "   p.x: " << p.x << "   p.y: " << p.y << std::endl;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	
  for( LandmarkObs& o: observations )
  {
    float min_dist = 99999999.9;      
    for( LandmarkObs& p: predicted)
    {
      float distance = dist(o.x, o.y, p.x, p.y); 
      if( distance < min_dist )
      {
        min_dist = distance;
        o.id = p.id;
        //std::cout << "min_dist: " << min_dist << std::endl;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) 
{
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
    // Standard deviations for x, y
   const double sig_x = std_landmark[0];
   const double sig_y = std_landmark[0];
  
  for( const LandmarkObs& o: observations )
  {
  	std::cout << "observations.x: " << o.x << "     observations.y: " << o.y <<std::endl;
  }
    
  // For each particle calculate the observation it makes
  // Transform landmarks to local measurements
  /*
  for(Particle &p: particles)
  {	
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();
   
    for(const Map::single_landmark_s& m_lm: map_landmarks.landmark_list)
    {
      double local_x = m_lm.x_f - p.x;
      double local_y = m_lm.y_f - p.y;      
      
      if(local_x < sensor_range && local_y < sensor_range)
      {
        // map to local
        double sensed_x = (cos(p.theta) * local_x) + (sin(p.theta) * local_y);
        double sensed_y = -(sin(p.theta) * local_x) + (cos(p.theta) * local_y);
      	
        std::cout << "p.id:" << p.id << "  p.theta: " << p.theta <<"   sensed_x: " << sensed_x  << "  sensed_y: " << sensed_y << "  m_lm.id_i: " << m_lm.id_i << std::endl;
      	p.associations.push_back(m_lm.id_i);
      	p.sense_x.push_back(sensed_x); //distance or position of landmark? 
      	p.sense_y.push_back(sensed_y);
      }
    }
  }
  exit(-1);
  */
  
  for(Particle &p: particles)
  {
    p.associations.clear();
    p.sense_x.clear();
    p.sense_y.clear();
    
    for(const Map::single_landmark_s& m_lm: map_landmarks.landmark_list)
    {
      if(dist(p.x, p.y, m_lm.x_f, m_lm.y_f) < sensor_range)
      {
        p.associations.push_back(m_lm.id_i);
        p.sense_x.push_back(m_lm.x_f); //distance or position of landmark? 
        p.sense_y.push_back(m_lm.y_f);
      }
    }
  }
  
  // create array of predictions
  for(Particle& p: particles)
  {	    
    vector<LandmarkObs> particle_prediction;
    for(const Map::single_landmark_s& m_lm: map_landmarks.landmark_list)
    {
      double local_x = m_lm.x_f - p.x;
      double local_y = m_lm.y_f - p.y;
      
      if(local_x < sensor_range && local_y < sensor_range)
      {
        // map to local
        double sensed_x = (cos(p.theta) * local_x) + (sin(p.theta) * local_y);
        double sensed_y = -(sin(p.theta) * local_x) + (cos(p.theta) * local_y);
        LandmarkObs lm{m_lm.id_i, sensed_x, sensed_y};
      	particle_prediction.push_back(lm);
      }
    }
    
    vector<LandmarkObs> observation_asso(observations);
   	dataAssociation(particle_prediction, observation_asso);
    
    // calculate normalization term
    double gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
    p.weight = 1;
    
    for( const LandmarkObs& o: observation_asso)
    {
      // calculate exponent
      int association_id = -1;
      for( unsigned int i = 0; i < particle_prediction.size(); ++i )
      {
        if( o.id == particle_prediction[i].id )
        {
          association_id = i;
          break;
        }        
      }
      double exponent;
      double delta_x = o.x - particle_prediction[association_id].x;
      double delta_y = o.y - particle_prediction[association_id].y;
      //std::cout << "delta_x: " << delta_x << "  delta_y: " << delta_y << std::endl;
      exponent = (pow(delta_x, 2) / (2 * pow(sig_x, 2))) + (pow(delta_y, 2) / (2 * pow(sig_y, 2)));
      
      // calculate weight using normalization terms and exponent
      double weight;
      weight = gauss_norm * exp(-exponent);
      p.weight *= weight; 
      //std::cout << p.id << ": exponent: " << exponent << " weight: " << weight << std::endl;
    }
     std::cout << p.id << ": weight: " << p.weight << std::endl;
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  double weight_sum = 0.0;
  for(Particle& p: particles)
  {
  	weight_sum += p.weight;
  }
  
  double max_normalized = 0.0;
  for(Particle& p: particles)
  {
  	p.weight = p.weight/weight_sum;
    if( p.weight > max_normalized) max_normalized = p.weight; 
  }
  
  
  double beta = 0.0;
  int N = particles.size();
  std::default_random_engine gen;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::uniform_int_distribution<int> dist_int(0.0, N-1);
  std::vector<Particle> particles_old(particles);
  
  for(Particle& p: particles)
  {
  	//std::cout <<  "resample  p.id: " << p.id << std::endl;
  }
  
  particles.clear();
  
  for( int i =0; i < N; ++i )
  {
    int index = dist_int(gen);
    beta = beta + dist(gen) * 2*max_normalized;
        
    while( particles[index].weight < beta )
    {
      beta = beta-particles[index].weight;
      index = (index+1) % N;
    }
    particles.push_back(particles_old[index]);		//weight
  }
  
  for(Particle& p: particles)
  {
  	std::cout << "resample  p.id: " << p.id << std::endl;
  }
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}