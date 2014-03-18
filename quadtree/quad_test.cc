#include <boost/test/unit_test.hpp>
#include <usml/types/types.h>
#include <iostream>

BOOST_AUTO_TEST_SUITE(quadtree_test)

using namespace std;
using namespace boost::unit_test;
using namespace usml::types;

BOOST_AUTO_TEST_CASE(quad_test){

  cout << "=== quadtree: build_test ===" << endl;
  quadtree* tree = new quadtree(0.0,USML_TEST_DIR "/types/test/wedge_cartesian.txt");
  tree->populate_tree();
  cout << "Tree built" << endl;
  BOOST_CHECK_EQUAL(tree->numberofleaves,8001);

  cout << "=== quadtree: traverse_test ===" << endl;
  tree->traverse(0,0,99,99);
  cout << "Tree traversed" << endl;
  BOOST_CHECK_EQUAL((tree->Search_Results.size())/3,199);

  cout << "=== quadtree: transform_test ===" << endl;
  std::vector<double> coords = tree->forward_transform(890.0,10.0,10.0,-200);
  coords = tree->inverse_transform(coords.at(1),coords.at(2),coords.at(0),-200);
  cout << "Coordinates transformed" << endl;
  BOOST_CHECK_CLOSE(coords.at(0),890.0,1e-3);
  BOOST_CHECK_CLOSE(coords.at(1),10.0,1e-3);
  BOOST_CHECK_CLOSE(coords.at(2),10.0,1e-3);
  coords.erase(coords.begin(),coords.end());

  delete tree;
}

BOOST_AUTO_TEST_SUITE_END()
