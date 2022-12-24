#include "remove_outliers.hpp"


void removeOutliers (Matcher::p_match p_matched[POINT_L], int32_t &p_matched_cnt){
    // do we have enough points for outlier removal?
    if (p_matched_cnt<=3)
        return;

    // input/output structure for triangulation
    delaunator::Delaunator d;
    int tmp = p_matched_cnt;
  
  
    // create copy of p_matched, init vector with number of support points
    // and fill triangle point vector for delaunay triangulation
    Vector<Matcher::p_match, POINT_L> p_matched_copy;  

    int32_t num_support[POINT_L][3];
//#pragma HLS ARRAY_PARTITION variable= num_support type=complete dim=2

    // initialize_loop: 
    for (int it = 0; it < p_matched_cnt; it++) {
        p_matched_copy.push_back(p_matched[it]);
        for(int j = 0; j < 3; j++)
        	num_support[it][j]=0;
        d.read_point(p_matched[it].u1c, p_matched[it].v1c);
    }
 
    // do triangulation (z=zero-based, n=neighbors, Q=quiet, B=no boundary markers)
    // attention: not using the B switch or using the n switch creates a memory leak (=> use valgrind!)
    
    d.delaunat();

    int32_t outlier_flow_tolerance = 5;
 
    // for all triangles do
    // accumulate_support: 
    for (int32_t i=0; i<d.triangles_cnt / 3; i++) {
        

        // extract triangle corner points
        int32_t p1 = d.triangles[i * 3 + 0];
        Matcher::p_match p1_match = p_matched_copy.fetch(p1);


        int32_t p2 = d.triangles[i * 3 + 1];
        Matcher::p_match p2_match = p_matched_copy.fetch(p2);      

        int32_t p3 = d.triangles[i * 3 + 2];
        Matcher::p_match p3_match = p_matched_copy.fetch(p3);

        // 1. corner disparity and flow
        float p1_flow_u = p1_match.u1c-p1_match.u1p;
        float p1_flow_v = p1_match.v1c-p1_match.v1p;

        // 2. corner disparity and flow
        float p2_flow_u = p2_match.u1c-p2_match.u1p;
        float p2_flow_v = p2_match.v1c-p2_match.v1p;

        // 3. corner disparity and flow
        float p3_flow_u = p3_match.u1c-p3_match.u1p;
        float p3_flow_v = p3_match.v1c-p3_match.v1p;


        bool flag1 = fabs(p1_flow_u-p2_flow_u)+fabs(p1_flow_v-p2_flow_v)<outlier_flow_tolerance; // consistency of 1. edge
        bool flag2 = fabs(p2_flow_u-p3_flow_u)+fabs(p2_flow_v-p3_flow_v)<outlier_flow_tolerance; // consistency of 2. edge
        bool flag3 = fabs(p1_flow_u-p3_flow_u)+fabs(p1_flow_v-p3_flow_v)<outlier_flow_tolerance; // consistency of 3. edge

        if(flag1 && flag3)
        	num_support[p1][0] += 2;
        else if(flag1 || flag3)
        	num_support[p1][0] += 1;

        if(flag1 && flag2)
        	num_support[p2][1] += 2;
        else if(flag1 || flag2)
        	num_support[p2][1] += 1;

        if(flag2 && flag3)
        	num_support[p3][2] += 2;
        else if(flag2 || flag3)
        	num_support[p3][2] += 1;
    } 

    // refill p_matched
    p_matched_cnt = 0;
    // refill: 
    for (int i=0; i<tmp; i++)
        if (num_support[i][0] + num_support[i][1] + num_support[i][2] >=4)
            p_matched[p_matched_cnt++] = p_matched_copy.fetch(i);
  

    
}
