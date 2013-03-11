// Created 12-Aug-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// g++ -o catalogAnalysis catalogAnalysis.cc quartic.cc -lboost_program_options -lgsl -llikely

/* = TODO =

- Implement option to size each ellipse to a specified surface brightness contour.
- Include bulge components from catalog (how?)
- Extend catalog to fainter objects.
- Compare shear correlation function with perfect deblending to other scenarios
- Remove likely dependency (using boost::random directly instead)
- Add fast rotated bounding box overlap test before full quartic test?

*/

#include "boost/program_options.hpp"
#include "boost/foreach.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/connected_components.hpp"

#include "likely/likely.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>

namespace po = boost::program_options;

std::vector<double> solve_quartic(double c0, double c1, double c2, double c3, double c4);

double pi = 4*std::atan(1);

class Ellipse {
public:
    Ellipse(double x, double y, double rmajor, double rminor, double phi)
    : _x(x), _y(y), _rmajor(rmajor), _rminor(rminor), _phi(phi), _overlapFlag(false)
    {
        assert(_rmajor >= _rminor);
        // Calculate shear parameters.
        _emag = (rmajor-rminor)/(rminor+rmajor);
        _e1 = _emag*std::cos(2*phi);
        _e2 = _emag*std::sin(2*phi);
        _r = 2*rminor*rmajor/(rminor+rmajor);
        // Calculate bilinear parameters.
        double tmp = 1 + _emag*_emag;
        double xtmp = tmp - 2*_e1;
        double ytmp = tmp + 2*_e1;
        double denom = x*x*xtmp + y*y*ytmp - 4*_e2*x*y - _r*_r;
        _a00 = -xtmp/denom;
        _a01 = 2*_e2/denom;
        _a11 = -ytmp/denom;
        _b0 = 2*(x*xtmp - 2*_e2*y)/denom;
        _b1 = 2*(-2*_e2*x + (1+_e1)*(1+_e1)*y + _e2*_e2*y)/denom;
    }
    bool isOverlapping(const Ellipse &other) const {
        
        // Check if bounding circles are non-overlapping.
        double dx = _x - other._x, dy = _y - other._y;
        double rsq = dx*dx + dy*dy;
        double rsum = _rmajor + other._rmajor;
        if(rsq > rsum*rsum) return false;
        
        // Quick circle test for overlap.
        rsum = _rminor + other._rminor;
        if(rsq < rsum*rsum) return true;
        
        // Use variables a00, etc, for the ellipse with the larger area so that there is only
        // one containment scenario that needs to be detected.
        double a00,a01,a11,b0,b1;
        double aa00,aa01,aa11,bb0,bb1;
        bool swap(_rmajor*_rminor < other._rmajor*other._rminor);
        if(swap) {
            aa00 = _a00; aa01 = _a01; aa11 = _a11; bb0 = _b0; bb1 = _b1;
            a00 = other._a00; a01 = other._a01; a11 = other._a11; b0 = other._b0; b1 = other._b1;
        }
        else {
            a00 = _a00; a01 = _a01; a11 = _a11; b0 = _b0; b1 = _b1;
            aa00 = other._a00; aa01 = other._a01; aa11 = other._a11; bb0 = other._b0; bb1 = other._b1;
        }

        // Precompute powers
        double a00sq(a00*a00),a01sq(a01*a01),a11sq(a11*a11),b0sq(b0*b0), b1sq(b1*b1);
        double aa00sq(aa00*aa00),aa01sq(aa01*aa01),aa11sq(aa11*aa11),bb0sq(bb0*bb0), bb1sq(bb1*bb1);
        double a01cu(a01*a01sq);
        
        double c0 = (-4*a01sq*a01sq + a11sq*aa00*b0sq + 2*a00*a11*b0*(aa01*b1 - a11*bb0) - 2*a01cu*(b1*bb0 + b0*bb1) + 
            a01sq*(8*a00*a11 + aa11*b0sq + 2*aa01*b0*b1 + aa00*b1sq + 2*a11*b0*bb0 + 2*a00*b1*bb1) + 
            a00sq*(-4*a11sq + aa11*b1sq - 2*a11*b1*bb1) - 
            2*a01*(a00*b1*(aa11*b0 + aa01*b1) + a11*(aa01*b0sq + aa00*b0*b1 - a00*(b1*bb0 + b0*bb1))))/4.;
             
        double c1 =  (a11*(-aa01sq + aa00*aa11)*b0sq - 2*a01cu*(4*aa01 + bb0*bb1) - a00sq*a11*(4*aa11 + bb1sq) + 
            a01sq*(4*a00*aa11 + 2*aa11*b0*bb0 - 2*aa01*b1*bb0 + a11*(4*aa00 + bb0sq) - 2*aa01*b0*bb1 + 2*aa00*b1*bb1 + 
            a00*bb1sq) + 2*a01*((aa01sq - aa00*aa11)*b0*b1 + a00*a11*(4*aa01 + bb0*bb1)) - 
            a00*((aa01sq - aa00*aa11)*b1sq + a11sq*(4*aa00 + bb0sq) + 
            2*a11*(aa11*b0*bb0 + aa00*b1*bb1 - aa01*(b1*bb0 + b0*bb1))))/2.;

        double c2 =   (-4*a00sq*aa11sq - aa01sq*aa11*b0sq + aa00*aa11sq*b0sq + 2*aa01*aa01sq*b0*b1 - 
            2*aa00*aa01*aa11*b0*b1 - aa00*aa01sq*b1sq + aa00sq*aa11*b1sq - 2*a00*aa11sq*b0*bb0 + 
            2*a00*aa01*aa11*b1*bb0 - a11sq*aa00*(4*aa00 + bb0sq) + 2*a00*aa01*aa11*b0*bb1 - 2*a00*aa01sq*b1*bb1 - 
            a00sq*aa11*bb1sq + a01sq*(-24*aa01sq + 8*aa00*aa11 + 3*aa11*bb0sq - 10*aa01*bb0*bb1 + 
            3*aa00*bb1sq) - 2*a01*(aa01sq*(b1*bb0 + b0*bb1) + aa00*aa11*(b1*bb0 + b0*bb1) - 
            2*aa01*(aa11*b0*bb0 + aa00*b1*bb1) - a00*(aa11*bb0*bb1 + aa01*(8*aa11 + bb1sq))) + 
            2*a11*(-((aa01*b0 - aa00*b1)*(aa01*bb0 - aa00*bb1)) + a01*(8*aa00*aa01 + aa01*bb0sq + aa00*bb0*bb1) + 
            a00*(4*aa01sq + 3*aa01*bb0*bb1 - 2*(4*aa00*aa11 + aa11*bb0sq + aa00*bb1sq))))/4.;

        double c3 = -((a11*aa00 - 2*a01*aa01 + a00*aa11)*(-4*aa01sq + aa11*bb0sq - 2*aa01*bb0*bb1 + aa00*(4*aa11 + bb1sq)))/2.;

        double c4 = -((aa01sq - aa00*aa11)*(4*aa01sq - aa11*bb0sq + 2*aa01*bb0*bb1 - aa00*(4*aa11 + bb1sq)))/4.;

        std::vector<double> roots = solve_quartic(c0,c1,c2,c3,c4);
        double lambdaMin(1),lambdaMax(1);
        BOOST_FOREACH(double t, roots) {
            // Calculate the (x,y) corresponding to this root of the quartic.
            double delta = -a01sq - 2*a01*aa01*t + a00*(a11 + aa11*t) + t*(a11*aa00 - aa01sq*t + aa00*aa11*t);
            if(0 == delta) continue;
            double x = (-((a11 + aa11*t)*(b0 + bb0*t)) + (a01 + aa01*t)*(b1 + bb1*t))/(2*delta);
            double y = ((a01 + aa01*t)*(b0 + bb0*t) - (a00 + aa00*t)*(b1 + bb1*t))/(2*delta);
            // Check that (x,y) is on the other ellipse
            if(std::fabs(swap ? getBilinear(x,y) : other.getBilinear(x,y)) > 1e-8*_r) continue;
            // Update the min,max values of lambda seen so far
            double lambda = swap ? other.getBilinear(x,y) : getBilinear(x,y);
            if(lambda < lambdaMin) lambdaMin = lambda;
            if(lambda > lambdaMax) lambdaMax = lambda;
            //std::cout << "(x,y,lam) = " << x << ',' << y << ',' << lambda << std::endl;
        }
        // lambdaMax < 0 means that we are completely contained with the other ellipse.
        return (lambdaMin <= 0 || lambdaMax < 0);
    }
    double getArea() const {
        return pi*_rmajor*_rminor;
    }
    double getBilinear(double x, double y) const {
        double result = -1 + _b0*x + _a00*x*x + _b1*y + 2*_a01*x*y + _a11*y*y;
        return (_a00 < 0) ? -result : +result;
    }
    double getFieldRadius() const {
        return std::sqrt(_x*_x + _y*_y);
    }
    void setOverlapFlag(bool value = true) {
        _overlapFlag = value;
    }
    bool getOverlapFlag() const {
        return _overlapFlag;
    }
    void dumpMath(std::ostream &os) const {
        os << "drawEllipse[" << _x << ',' << _y << ',' << _rmajor << ',' << _rminor << ',' << _phi
            << ',' << (getOverlapFlag() ? "True":"False") << "]" << std::endl;
    }
    void dumpRoot(std::ostream &os) const {
        os << _x << ' ' << _y << ' ' << _rmajor << ' ' << _rminor << ' ' << _emag
            << ' ' << (_overlapFlag ? 1:0) << std::endl;
    }
private:
    // Ellipse center
    double _x,_y;
    // Semi-major/minor axes and anti-cw rotation of semi-major axis from x axis in radians.
    double _rmajor,_rminor,_phi;
    // Shear parameters
    double _emag, _e1, _e2, _r;
    // Bilinear parameters (encoding position + shape)
    double _a00, _a01, _a11, _b0, _b1;
    // Has an overlap been detected for this ellipse?
    bool _overlapFlag;
};

int main(int argc, char *argv[]) {
    // Configure cmd-line options
    po::options_description cli("Geometrical overlaps calculator");
    int seed;
    double fieldCenterRA,fieldCenterDec,fieldRadius,density,galaxyRadius,padFraction,psfSize;
    std::string mathOutput,rootOutput,graphOutput,catalogName;
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("field-center-ra", po::value<double>(&fieldCenterRA)->default_value(0.5),
            "RA of field center in deg.")
        ("field-center-dec", po::value<double>(&fieldCenterDec)->default_value(-0.5),
            "Dec of field center in deg.")
        ("field-radius", po::value<double>(&fieldRadius)->default_value(1),
            "Size of field in arcminutes.")
        ("pad-fraction", po::value<double>(&padFraction)->default_value(0.1),
            "Fractional padding to remove edge effects.")
        ("density", po::value<double>(&density)->default_value(100),
            "Galaxies per square arcminute (converted to catalog i-band cut if < 0).")
        ("galaxy-radius", po::value<double>(&galaxyRadius)->default_value(1),
            "Galaxy radius in arcseconds (scales half-light radius if < 0).")
        ("psf-size", po::value<double>(&psfSize)->default_value(0),
            "Amount to add in quadrature to semi-major,minor axis lengths in arcsecs.")
        ("math-output", po::value<std::string>(&mathOutput)->default_value(""),
            "Filename to save mathematica output.")
        ("root-output", po::value<std::string>(&rootOutput)->default_value(""),
            "Filename to save root tree data as plain text.")
        ("graph-output", po::value<std::string>(&graphOutput)->default_value(""),
            "Filename to save overlap graph analysis.")
        ("catalog-name", po::value<std::string>(&catalogName)->default_value(
            "/Volumes/Data/desc/catalogs/lsstSimCatalog.dat"),
            "Name of input catalog to read.")
        ("seed", po::value<int>(&seed)->default_value(123),
            "Random seed to use.")
        ("random-xy","Randomize (x,y) positions.")
        ("circular", "All galaxies are circular.")
        ;
    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose")),randomXY(vm.count("random-xy")),
        circular(vm.count("circular"));
    bool doGraphAnalysis = graphOutput.length() > 0;
    
    // Initialize random numbers.
    likely::RandomPtr random(likely::Random::instance());
    random->setSeed(seed);

    // Open catalog file and skip header line.
    std::ifstream cat(catalogName.c_str());
    std::string header;
    std::getline(cat,header);

    double galaxyArea = pi*galaxyRadius*galaxyRadius;
    double fiducialRadius = fieldRadius/(1 + padFraction);
    double psfSizeSq = psfSize*psfSize/3600.;
    
    // Create the galaxies in this field.
    std::vector<Ellipse> galaxies;
    if(density > 0) {
        // Generate galaxies at the requested density.
        int nGalaxy = std::floor(pi*fieldRadius*fieldRadius*density+0.5);
        for(int index = 0; index < nGalaxy; ++index) {
            double r = std::sqrt(random->getUniform())*fieldRadius;
            double theta = 2*pi*random->getUniform();
            double x = r*std::cos(theta);
            double y = r*std::sin(theta);
            double ratio = 1;
            double phi = 0;
            if(!circular) {
                // A ratio in the range [0.2,1] gives an ellipticity in the range [0,2/3]
                ratio = (1+4*random->getUniform())/5.;
                phi = random->getUniform()*2*pi;
            }
            // Calculate major,minor axes in arcmins from area and ratio = major/minor
            double rminor = std::sqrt(galaxyArea*ratio/pi)/60.;
            double rmajor = rminor/ratio;
            galaxies.push_back(Ellipse(x,y,rmajor,rminor,phi));
        }
    }
    else {
        // Estimate i-band cut for requested density.
        double icut = 1.40095*std::log(1.22248e6*(-density));
        // Loop over catalog entries.
        long id;
        int lines(0);
        // Per-galaxy columns (r_ab appears twice: here and in the per-filter section... to be fixed)
        double ra,dec,redshift,r_ab,absmag_r_total;
        // Per-component (bulge/disk) columns
        double BulgeHalfLightRadius,DiskHalfLightRadius,pa_bulge,pa_disk,magnorm_bulge,magnorm_disk,
            sedid_bulge,sedid_disk,a_b,a_d,b_b,b_d;
        std::string sedname_bulge,sedname_disk;
        // Per-filter (ugrizy) colunns
        double u_ab,g_ab,/*r_ab,*/i_ab,z_ab,y_ab;
        while(cat.good()) {
            cat >> id>>ra>>dec>>redshift>>r_ab>>absmag_r_total
                >> BulgeHalfLightRadius>>DiskHalfLightRadius>>pa_bulge>>pa_disk>>magnorm_bulge>>magnorm_disk
                >> sedid_bulge>>sedid_disk>>sedname_bulge>>sedname_disk>>a_b>>a_d>>b_b>>b_d
                >> u_ab>>g_ab>>r_ab>>i_ab>>z_ab>>y_ab;
            if(!cat.good()) break;
            lines++;
            // Does this pass our magnitude limit?
            if(i_ab > icut) continue;
            // Convert (ra,dec) in degrees to arcmins relative to the field center.
            double x = (ra - fieldCenterRA)*60;
            double y = (dec - fieldCenterDec)*60;
            // Is this object within our field.
            if(x*x + y*y > fieldRadius*fieldRadius) continue;
            // Create this galaxy.
            double ratio = 1;
            double phi = 0;
            if(!circular) {
                // Ignore galaxies with no disk component
                if(a_d == 0) continue;
                // The catalog is roughly flat in ratio = b_d/a_d
                ratio = b_d/a_d;
                phi = random->getUniform()*2*pi;
            }
            // Calculate the galaxy major and minor radii in arcmins.
            double rminor,rmajor;
            if(galaxyRadius < 0) {
                rminor = b_d/60.;
            }
            else {
                rminor = std::sqrt(galaxyArea*ratio/pi)/60.;
            }
            rmajor = rminor/ratio;
            if(psfSizeSq > 0) {
                rminor = std::sqrt(rminor*rminor + psfSizeSq);
                rmajor = std::sqrt(rmajor*rmajor + psfSizeSq);
            }
            if(galaxyRadius < 0) {
                rminor *= -galaxyRadius;
                rmajor *= -galaxyRadius;
            }
            // Add this galaxy to the field.
            galaxies.push_back(Ellipse(x,y,rmajor,rminor,phi));
        }
    }
    
    // Calculate overlaps and build adjacency graph.
    int nGalaxy = galaxies.size();
    boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph;
    for(int index1 = 0; index1 < nGalaxy-1; ++index1) {
        for(int index2 = index1+1; index2 < nGalaxy; ++index2) {
            // Is there any point in checking for an overlap?
            if(!doGraphAnalysis && galaxies[index1].getOverlapFlag() &&
                galaxies[index2].getOverlapFlag()) continue;
            if(galaxies[index1].isOverlapping(galaxies[index2])) {
                galaxies[index1].setOverlapFlag();
                galaxies[index2].setOverlapFlag();
                // Add this overlap to our undirected graph.
                if(doGraphAnalysis) boost::add_edge(index1,index2,graph);
            }
        }
    }
    
    // Analyze overlap graph.
    if(doGraphAnalysis) {
        // Identify clusters as connected subgraphs.
        std::vector<int> clusters(boost::num_vertices(graph));
        int nClusters = boost::connected_components(graph, &clusters[0]);
        // Count the number of galaxies in each cluster.
        std::vector<int> clusterSize(nClusters,0);
        BOOST_FOREACH(int cluster, clusters) clusterSize[cluster]++;
        // Sort clusters in increasing size.
        std::sort(clusterSize.begin(),clusterSize.end());
        // Dump the size of each cluster and the cummulative fraction of galaxies
        // in clusters at least this size.
        std::ofstream out(graphOutput.c_str());
        out << "csize frac" << std::endl;
        int ncum(0);
        BOOST_FOREACH(int size, clusterSize) {
            ncum += size;
            out << size << ' ' << (double)ncum/nGalaxy << std::endl;
        }
        out.close();
    }
    
    // Calculate overlap statistics.
    double ntotal,noverlap;
    for(int index = 0; index < nGalaxy; ++index) {
        if(galaxies[index].getFieldRadius() > fiducialRadius) continue;
        ntotal++;
        if(galaxies[index].getOverlapFlag()) noverlap++;
    }
    if(density < 0) density = ntotal/(pi*fiducialRadius*fiducialRadius);
    std::cout << density << ' ' << noverlap/ntotal << std::endl;
    
    // Dump the field to a mathematica file if requested.
    if(0 < mathOutput.length()) {
        std::ofstream out(mathOutput.c_str());
        out << "{ Black, Dashed, Circle[{0,0}," << fieldRadius
            << "], Dashing[{}], Circle[{0,0}," << fiducialRadius
            << "], Opacity[0.25], Blue," << std::endl;
        for(int index = 0; index < nGalaxy; ++index) {
            galaxies[index].dumpMath(out);
            if(index < nGalaxy-1) out << ',';
        }
        out << "}" << std::endl;
        out.close();
    }
    
    // Dump the field to a root output file if requested.
    if(0 < rootOutput.length()) {
        std::ofstream out(rootOutput.c_str());
        out << "x y a b e overlap" << std::endl;
        for(int index = 0; index < nGalaxy; ++index) {
            galaxies[index].dumpRoot(out);
        }
        out.close();
    }

    cat.close();
    return 0;
}
