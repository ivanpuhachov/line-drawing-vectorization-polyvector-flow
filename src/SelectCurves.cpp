#include "SelectCurves.h"

typedef std::set<std::pair<int, int>> PixelCoverage;

struct allowedCurve
{
    int id;
    Chain chain;
    std::set<std::pair<int, int>> curveCoverage;
    double curveLength;
    double curveEnergy;
};

MyPolyline flipXY(const MyPolyline& poly)
{
    MyPolyline result(poly.size());
    for (int i = 0; i < poly.size(); ++i)
        result[i] = Eigen::Vector2d(poly[i].y(), poly[i].x());
    return result;
}

std::set<std::pair<int, int>> computeCurvesCoverage(const std::vector<allowedCurve>& vectorAllowedCurves, std::set<int> selectionId) {
    std::set<std::pair<int, int>> totalCoverage;
    for (const allowedCurve& ac : vectorAllowedCurves) {
        if (selectionId.find(ac.id)!=selectionId.end()){
            totalCoverage.insert(ac.curveCoverage.begin(), ac.curveCoverage.end());
        }
    }
    return totalCoverage;
}

std::vector<Chain> selectAllCurves(
        const G& g,
        const std::vector<Chain>& currentChains,
        const cv::Mat& origMask,
        const FrameFieldFlow& fff
) {
    std::vector <allowedCurve> availableCurves;
    int id=0;
    for (const Chain& c : currentChains) {
        availableCurves.push_back(
                {id,
                 c,
                 computeCoverage(g, chainToVerticesSeq(c), origMask, false),
                 chainLength(c, g),
                 fff.fullEnergy(flipXY(chainToPolyline(c, g)))});
        id++;
    }
    std::cout << "\n\n\nALL CURVES: \n";
    for (id=0; id<availableCurves.size(); id++){
        std::cout << "Curve " << availableCurves[id].id << ", length: " << availableCurves[id].curveLength << ", energy: " << availableCurves[id].curveEnergy;
        printChain(availableCurves[id].chain, g, false);
    }
    return currentChains;
}

std::vector<Chain> selectBestCurves(
        const G& g,
        const std::vector<Chain>& currentChains,
        const cv::Mat& origMask,
        const FrameFieldFlow& fff
        ) {
    std::vector <allowedCurve> availableCurves;
    int id=0;
    for (const Chain& c : currentChains) {
        availableCurves.push_back(
                {id,
                 c,
                 computeCoverage(g, chainToVerticesSeq(c), origMask, false),
                 chainLength(c, g),
                 fff.fullEnergy(flipXY(chainToPolyline(c, g)))});
        id++;
    }
    // sort in descending order
    std::sort(availableCurves.begin(), availableCurves.end(),
              [](const auto& i, const auto& j) { return i.curveEnergy > j.curveEnergy; } );
    for (id=0; id<availableCurves.size(); id++){
        std::cout << "Curve " << availableCurves[id].id << ", length: " << availableCurves[id].curveLength << ", energy: " << availableCurves[id].curveEnergy;
        printChain(availableCurves[id].chain, g, true);
    }
//    std::set<allowedCurve> selectedCurves(availableCurves.begin(), availableCurves.end());
    std::set<int> selectedCurvesIdsSet;
    for (id=0; id<availableCurves.size(); id++){
        selectedCurvesIdsSet.insert(id);
    }

    for (const allowedCurve& curve : availableCurves){
        std::cout <<"\n";
        std::cout << "Working with curve " << curve.id << "\n";
        double curveEnergy = fff.fullEnergy(flipXY(chainToPolyline(curve.chain, g)));
//        if (curve.id==5){
//            curveEnergy = fff.fullEnergy2(chainToPolyline(curve.chain, g));
//        }
        std::cout << "This curve energy is " << curveEnergy << "\n";
        selectedCurvesIdsSet.erase(curve.id);
        PixelCoverage othersCoverage = computeCurvesCoverage(availableCurves, selectedCurvesIdsSet);
        PixelCoverage residualCoverage;
        std::set_difference(curve.curveCoverage.begin(), curve.curveCoverage.end(), othersCoverage.begin(), othersCoverage.end(),
                            std::inserter(residualCoverage, residualCoverage.end()));
        std::cout << "This curve adds " << residualCoverage.size() << " new pixels. \n";
        if (residualCoverage.size()>2){
            std::cout << "It's bigger than threshold, return curve " << curve.id << " to the set of selected.\n";
            selectedCurvesIdsSet.insert(curve.id);
        }
    }

    std::vector<Chain> result;
    for (const allowedCurve& ac : availableCurves){
        if (selectedCurvesIdsSet.find(ac.id)!=selectedCurvesIdsSet.end()){
            result.push_back(ac.chain);
        }
    }

    return result;
}

double polylineLength(const MyPolyline& poly){
    double length = 0.0;
    for (size_t i=0; i<poly.size()-1; ++i){
        length+= (poly[i] - poly[i-1]).norm();
    }
    return length;
}

struct allowedPoly
{
    int id;
    MyPolyline poly;
    double polyLength;
    double polyEnergy;
    std::pair<int, int> labeledPointEnds;
};

std::vector<MyPolyline> selectBestPolylines(
        const G& g,
        const std::vector<MyPolyline>& currentPolys,
        const std::vector<MyPolyline>& initialPolys,
        const std::vector<std::vector<vertex_descriptor>>& chains,
        const FrameFieldFlow& fff,
        std::vector<std::pair<int, int>> allowedChainToLabeledPointPair,
        const std::vector<int>& pts_types,
        std::map<std::pair<int,int>, std::vector<size_t>> pixelToChainCoveringMap,
        std::vector<bool>& selection
) {
    selection.clear();
    std::vector <allowedPoly> availablePolys;
    int id=0;
    for (size_t i_poly=0; i_poly<currentPolys.size(); ++i_poly) {
        MyPolyline poly = currentPolys[i_poly];
//        double poly_energy = fff.fullEnergy(flipXY(poly));
        double poly_energy = fff.fullEnergy_white(flipXY(poly), flipXY(initialPolys[i_poly]), chains[i_poly]);
        double poly_length = polylineLength(poly);
        std::pair<int, int> ends = allowedChainToLabeledPointPair[i_poly];
        availablePolys.push_back({
            id,
            poly,
            poly_length,
            poly_energy,
            ends
            });
        id++;
    }
    size_t n_labeledpoints = pts_types.size();
    size_t n_polys = currentPolys.size();

    for (size_t i=0; i<n_labeledpoints; ++i){
        std::cout << "Keypoint: " << i << " has type " << pts_types[i] << "\n";
    }

    for (id=0; id<availablePolys.size(); id++){
        std::cout << "Poly: " << availablePolys[id].id << "; with len " <<
        availablePolys[id].polyLength << "; with energy " << availablePolys[id].polyEnergy << "; it connects "
        << availablePolys[id].labeledPointEnds.first <<" and " << availablePolys[id].labeledPointEnds.second << "\n";
    }

    std::unique_ptr<GRBEnv> env_ = std::unique_ptr<GRBEnv>(new GRBEnv());
    GRBModel model(*env_);

    // creating variables
    auto vars = std::unique_ptr<GRBVar[]>(model.addVars(n_polys, GRB_BINARY));
    model.update();

// valence constraints
    for (size_t i_labeled=0; i_labeled<n_labeledpoints; ++i_labeled){
        GRBLinExpr lhs = 0;
        int maxvalence = 0;
        for (size_t i_poly=0; i_poly<n_polys; ++i_poly){
            std::pair<int, int> poly_ends = allowedChainToLabeledPointPair[i_poly];
            if ((poly_ends.first==i_labeled)||(poly_ends.second==i_labeled)){
                lhs += vars[i_poly];
                maxvalence++;
            }

        }
        int rhs = 0;
        if (pts_types[i_labeled]==endpoint){
            rhs = std::min(maxvalence, 1);
            if (maxvalence<1){
                std::cout << "POINT " << i_labeled << " DOES NOT PASS type: " << pts_types[i_labeled] << ", maxvalence: " << maxvalence << "\n";
            }
        }
        if (pts_types[i_labeled]==sharpcorner){
            rhs = std::min(maxvalence, 2);
            if (maxvalence<2){
                std::cout << "POINT " << i_labeled << " DOES NOT PASS type: " << pts_types[i_labeled] << ", maxvalence: " << maxvalence << "\n";
            }
        }
        if (pts_types[i_labeled]==junction){
            rhs = std::min(maxvalence, 3);
            if (maxvalence<3){
                std::cout << "POINT " << i_labeled << " DOES NOT PASS type: " << pts_types[i_labeled] << ", maxvalence: " << maxvalence << "\n";
            }
        }
        model.addConstr(lhs, GRB_GREATER_EQUAL, rhs);
    }

// coverage constraints
    for ( const auto &keyvalupair : pixelToChainCoveringMap) {
        GRBLinExpr lhs = 0;
        for (const size_t i_poly : keyvalupair.second){
            lhs += vars[i_poly];
        }
        model.addConstr(lhs, GRB_GREATER_EQUAL, 1);
    }
//    model.addConstr(vars[3] + vars[1] + vars[5], GRB_GREATER_EQUAL, 2);
//    model.addConstr(vars[3] + vars[1] + vars[5], GRB_GREATER_EQUAL, 2);
//    model.addConstr(vars[2] + vars[4] + vars[6], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[3] + vars[4], GRB_GREATER_EQUAL, 1);


// coverage constraints
//    model.addConstr(vars[3] + vars[4], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[6] + vars[4], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[0], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[1], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[2], GRB_GREATER_EQUAL, 1);
//    model.addConstr(vars[0] + vars[4], GRB_GREATER_EQUAL, 1);

    GRBLinExpr loss;
    for (int i=0; i<currentPolys.size(); ++i){
        loss += vars[i] * availablePolys[i].polyEnergy;
    }
    model.setObjective(loss, GRB_MINIMIZE);

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);

    if (optimstatus == GRB_OPTIMAL
        || optimstatus== GRB_SUBOPTIMAL) {

        std::vector<MyPolyline> solvedPolys;
        std::vector<double> solution(currentPolys.size(), 0.0);
        for (int i = 0; i < currentPolys.size(); i++) {
            solution[i] = vars[i].get(GRB_DoubleAttr_X);
            std::cout << "Solution " << i << " = " << solution[i] << "\n";
            if (solution[i]==1) {
                solvedPolys.push_back(currentPolys[i]);
                selection.push_back(true);
            }
            else{
                selection.push_back(false);
            }
        }
        return solvedPolys;
    } else {
        return currentPolys;
    }
    return currentPolys;


//    // sort in descending order
//    std::sort(availableCurves.begin(), availableCurves.end(),
//              [](const auto& i, const auto& j) { return i.curveEnergy > j.curveEnergy; } );
//    for (id=0; id<availableCurves.size(); id++){
//        std::cout << "Curve " << availableCurves[id].id << ", length: " << availableCurves[id].curveLength << ", energy: " << availableCurves[id].curveEnergy;
//        printChain(availableCurves[id].chain, g, true);
//    }
////    std::set<allowedCurve> selectedCurves(availableCurves.begin(), availableCurves.end());
//    std::set<int> selectedCurvesIdsSet;
//    for (id=0; id<availableCurves.size(); id++){
//        selectedCurvesIdsSet.insert(id);
//    }
//
//    for (const allowedCurve& curve : availableCurves){
//        std::cout <<"\n";
//        std::cout << "Working with curve " << curve.id << "\n";
//        double curveEnergy = fff.fullEnergy(chainToPolyline(curve.chain, g));
////        if (curve.id==5){
////            curveEnergy = fff.fullEnergy2(chainToPolyline(curve.chain, g));
////        }
//        std::cout << "This curve energy is " << curveEnergy << "\n";
//        selectedCurvesIdsSet.erase(curve.id);
//        PixelCoverage othersCoverage = computeCurvesCoverage(availableCurves, selectedCurvesIdsSet);
//        PixelCoverage residualCoverage;
//        std::set_difference(curve.curveCoverage.begin(), curve.curveCoverage.end(), othersCoverage.begin(), othersCoverage.end(),
//                            std::inserter(residualCoverage, residualCoverage.end()));
//        std::cout << "This curve adds " << residualCoverage.size() << " new pixels. \n";
//        if (residualCoverage.size()>2){
//            std::cout << "It's bigger than threshold, return curve " << curve.id << " to the set of selected.\n";
//            selectedCurvesIdsSet.insert(curve.id);
//        }
//    }
//
//    std::vector<Chain> result;
//    for (const allowedCurve& ac : availableCurves){
//        if (selectedCurvesIdsSet.find(ac.id)!=selectedCurvesIdsSet.end()){
//            result.push_back(ac.chain);
//        }
//    }
//
//    std::cout << "\n\n\nALL CURVES: \n";
//    for (id=0; id<availableCurves.size(); id++){
//        std::cout << "Curve " << availableCurves[id].id << ", length: " << availableCurves[id].curveLength << ", energy: " << availableCurves[id].curveEnergy;
//        printChain(availableCurves[id].chain, g, false);
//    }
//    std::cout << "\n\n\n";
//    return result;
}
