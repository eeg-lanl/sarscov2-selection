#ifndef GENEALOGY_HPP_
#define GENEALOGY_HPP_

#include <map>
#include <list>

#include "particle.hpp"
#include "aux.hpp"

typedef std::list<Particle> ParticlePath;

/** take a path of Particles, extract their stored States
 * and concatenate them together in a ParticlePath
 */
StatePath makeStatePath(const ParticlePath & ppath);

class Genealogy : public Printable {
public:
  /** add a new layer to the Genealogy.
   * Extinct lineages are removed.
   * @param particles the current particle swarm
   * @param counts the number of times a particle is re-sampled
   */
  void update(const std::vector<Particle*> & particles,
      const std::map<int, int> & counts);
  /** Check if the Genealogy is empty.
   * This checks both if particle_idx_maps is empty and
   * if ancestries is empty
   */
  bool empty() const; // test if the genealogy is empty
  /** Clear the content of the Genealogy.
   */
  void clear();
  /** Sample a path from the Genealogy.
   */
  ParticlePath sampleParticlePath(Rng & rng) const;
  StatePath sampleStatePath(Rng & rng) const;
  /** Print a summary of the Genealogy.
   * Currently this is for testing...
   */
  void print(std::ostream & os) const;
protected:
  typedef std::map<int, Particle> ParticleIdxMap;
  typedef std::map<int, int> AncestorIdxMap;
  std::list<ParticleIdxMap> particle_idx_maps;
  std::list<AncestorIdxMap> ancestries;
};



#endif
