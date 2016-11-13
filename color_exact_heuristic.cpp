/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 * Usa ideias e código de Rafael Arakaki e Flávio Keidi Miyazawa
 ******************************************************************************/

/* IMPLEMENTE AS FUNCOES INDICADAS
 * DIGITE SEU RA: 147775
 * SUBMETA SOMENTE ESTE ARQUIVO */

#include <iostream>
#include <gurobi_c++.h>
#include <float.h>
#include <stack>
#include <deque>
#include <vector>
#include <map>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "color_exact_heuristic.h"

bool vcomp(GRBVar, GRBVar) throw(GRBException);
int colorNaive(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit);

//------------------------------------------------------------------------------
int colorExact(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* SUBSTITUA O CONTEÚDO DESTA FUNÇÃO POR SUA IMPLEMENTAÇÃO DO ALGORITMO EXATO.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DA FUNÇÃO. */
{
  char name[1000];
  int nodeCount = 0;
  int status;
  for (NodeIt v(gd.g); v!=INVALID; ++v) nodeCount++;
  vector<GRBVar> y(nodeCount);
  ListGraph::NodeMap<vector<GRBVar>> x(gd.g);

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_MIPFocus, 1);
  upperBound = nodeCount;
  env.set(GRB_DoubleParam_Cutoff, upperBound);
  model.set(GRB_StringAttr_ModelName, "k Coloring with GUROBI"); // name to the problem
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

  // add y vars
  for (int i = 0; i < nodeCount; i++) {
    sprintf(name,"y_%d",i);
    y[i] = model.addVar(0.0, 1.0, 1.0 ,GRB_BINARY,name);
  }

  // add x vars (node x, color i)
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    x[v] = vector<GRBVar>(nodeCount);
    for(int i = 0; i < nodeCount; i++){
      sprintf(name,"x_%s,%d",gd.vname[v].c_str(), i);
      x[v][i] = model.addVar(0.0, 1.0, 0.0 ,GRB_BINARY, name);
    }
  }

  // the sum for Xik, for each i, vertex, must be equal to 1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    GRBLinExpr expr;
    for (int i =0; i<nodeCount;i++) expr += x[v][i];
    model.addConstr(expr == 1 );
  }

  // only set vertex v to use color i, if color i is really =1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    for (int i =0; i<nodeCount;i++){
      GRBLinExpr expr;
      expr = x[v][i] - y[i];
      model.addConstr(expr <= 0 );
    }
  }

  // only set vertex v to use color i, if color i is really =1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    for (int i =0; i<nodeCount;i++){
      GRBLinExpr expr;
      expr = x[v][i] - y[i];
      model.addConstr(expr <= 0 );
    }
  }
  // each adjacent vertices can't have the same color
  for (EdgeIt e(gd.g); e!=INVALID; ++e) {
    for(int i = 0; i < nodeCount; i++){
      GRBLinExpr expr;
      expr = x[gd.g.u(e)][i] + x[gd.g.v(e)][i];
      model.addConstr(expr <= 1 );
    }
  }
  try {
    if (timeLimit >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,timeLimit);
    model.update(); // run update to use model inserted variables
    model.optimize();

    map <int, int> colorMap;
    int j = 1;
    for(int i = 0; i < nodeCount; i++){
      if(BinaryIsOne(y[i].get(GRB_DoubleAttr_X))){
        colorMap[i] = j;
        j++;
      }
    }
    for (NodeIt v(gd.g); v!=INVALID; ++v) {
      for (int i =0; i<nodeCount;i++){
        if(BinaryIsOne(x[v][i].get(GRB_DoubleAttr_X))){
          color[v] = colorMap[i];
        }
      }
    }
    lowerBound = model.get(GRB_DoubleAttr_ObjBound);
    upperBound = model.get(GRB_DoubleAttr_ObjVal);
    status = model.get(GRB_IntAttr_Status);
  } catch (...){

  }
  // 2 = OPTIMAL
  return status == 2 ? 1 : 0;
}
//------------------------------------------------------------------------------
int colorHeuristic(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* SUBSTITUA O CONTEÚDO DESTA FUNÇÃO POR SUA IMPLEMENTAÇÃO DA HEURISTICA.
 * ENTRETANTO, NÃO ALTERE A ASSINATURA DA FUNÇÃO. */
{
  int remainingTime = timeLimit;
  const clock_t begin_time = clock();
  char name[1000];
  int nodeCount = 0;
  int status;
  for (NodeIt v(gd.g); v!=INVALID; ++v) nodeCount++;
  vector<GRBVar> y(nodeCount);
  ListGraph::NodeMap<vector<GRBVar>> x(gd.g);

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntParam_MIPFocus, 1);
  upperBound = nodeCount;
  env.set(GRB_DoubleParam_Cutoff, upperBound);
  model.set(GRB_StringAttr_ModelName, "k Coloring with GUROBI"); // name to the problem
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

  // add y vars
  for (int i = 0; i < nodeCount; i++) {
    sprintf(name,"y_%d",i);
      y[i] = model.addVar(0.0, 1.0, 1.0 ,GRB_CONTINUOUS,name);
  }

  // add x vars (node x, color i)
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    x[v] = vector<GRBVar>(nodeCount);
    for(int i = 0; i < nodeCount; i++){
      sprintf(name,"x_%s,%d",gd.vname[v].c_str(), i);
      x[v][i] = model.addVar(0.0, 1.0, 0.0 ,GRB_CONTINUOUS, name);
    }
  }

  // the sum for Xik, for each i, vertex, must be equal to 1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    GRBLinExpr expr;
    for (int i =0; i<nodeCount;i++) expr += x[v][i];
    model.addConstr(expr == 1 );
  }

  // only set vertex v to use color i, if color i is really =1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    for (int i =0; i<nodeCount;i++){
      GRBLinExpr expr;
      expr = x[v][i] - y[i];
      model.addConstr(expr <= 0 );
    }
  }

  // only set vertex v to use color i, if color i is really =1
  for (NodeIt v(gd.g); v!=INVALID; ++v) {
    for (int i =0; i<nodeCount;i++){
      GRBLinExpr expr;
      expr = x[v][i] - y[i];
      model.addConstr(expr <= 0 );
    }
  }
  // each adjacent vertices can't have the same color
  for (EdgeIt e(gd.g); e!=INVALID; ++e) {
    for(int i = 0; i < nodeCount; i++){
      GRBLinExpr expr;
      expr = x[gd.g.u(e)][i] + x[gd.g.v(e)][i];
      model.addConstr(expr <= 1 );
    }
  }
  try {
    model.update();
    if (timeLimit >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,remainingTime);
    // Solve the RELAXED PROBLEM
    model.optimize();
    lowerBound = model.get(GRB_DoubleAttr_ObjVal);
    bool shouldTry = true;
    for (int iter = 0; iter < 1000 && shouldTry; ++iter)
    {
        while(true){
            // create a list of fractional variables, sorted in order of
            // increasing distance from the relaxation solution to the nearest
            // integer value
            deque<GRBVar> fractional;
            for (NodeIt v(gd.g); v!=INVALID; ++v) {
                for (size_t j = 0; j < y.size(); j++)
                {
                    double sol = fabs(x[v][j].get(GRB_DoubleAttr_X));
                    if (fabs(sol - floor(sol + 0.5)) > 1e-5)
                    {
                        fractional.push_back(x[v][j]);
                    }
                }
                
            }
            
            cout << "Iteration " << iter << ", obj " <<
            model.get(GRB_DoubleAttr_ObjVal) << ", fractional " <<
            fractional.size() << endl;
            
            if (fractional.size() == 0)
            {
                shouldTry = false;
                cout << "Found feasible solution - objective " <<
                model.get(GRB_DoubleAttr_ObjVal) << endl;
                break;
            }
            
            int randomIndex = rand() % fractional.size();
            GRBVar v = fractional[randomIndex];
            v.set(GRB_DoubleAttr_LB, 1.0);
            v.set(GRB_DoubleAttr_UB, 1.0);
            
            model.optimize();
            
            // Check optimization result
            // if relaxation is infeasible
            if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
            {
                // Reset lower and upper bound for each variable
                for (NodeIt v(gd.g); v!=INVALID; ++v) {
                    for (size_t j = 0; j < y.size(); j++)
                    {
                        GRBVar var = x[v][j];
                        var.set(GRB_DoubleAttr_LB, 0.0);
                        var.set(GRB_DoubleAttr_UB, 1.0);
                    }
                    
                }
                cout << "Relaxation is infeasible" << endl;
                break;
            }
        }
        
    }

    map <int, int> colorMap;
    int j = 1;
    for(int i = 0; i < nodeCount; i++){
      if(BinaryIsOne(y[i].get(GRB_DoubleAttr_X))){
        colorMap[i] = j;
        j++;
      }
    }
    for (NodeIt v(gd.g); v!=INVALID; ++v) {
      for (int i =0; i<nodeCount;i++){
//          cout << "x_" << gd.vname[v] << "_" << i << ":" << x[v][i].get(GRB_DoubleAttr_X) << endl;
        if(x[v][i].get(GRB_DoubleAttr_X)>0){
          color[v] = colorMap[i];
        }
      }
    }
      upperBound = model.get(GRB_DoubleAttr_ObjVal);
  } catch(GRBException e){
      cout << e.getMessage() << endl;
  }
  return 0;
}
//------------------------------------------------------------------------------
int colorNaive(GraphData& gd, NodeIntMap& color, int& lowerBound, int& upperBound, int timeLimit)
/* Algoritmo ingênuo para o problema de coloração de vértices */
{
   lowerBound = 1;
   int next = lowerBound;

	for(NodeIt i(gd.g); i != INVALID; ++i){
      color[i] = next;
      next++;
	}
	next--;
	upperBound = next;

	return 0;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
bool vcomp(GRBVar v1,
           GRBVar v2) throw(GRBException)
{
  double sol1 = fabs(v1.get(GRB_DoubleAttr_X));
  double sol2 = fabs(v2.get(GRB_DoubleAttr_X));
  double frac1 = fabs(sol1 - floor(sol1 + 0.5));
  double frac2 = fabs(sol2 - floor(sol2 + 0.5));
  return (frac1 < frac2);
}
//------------------------------------------------------------------------------
