#= 
* This file is part of OpenModelica.
*
* Copyright (c) 1998-CurrentYear, Open Source Modelica Consortium (OSMC),
* c/o Linköpings universitet, Department of Computer and Information Science,
* SE-58183 Linköping, Sweden.
*
* All rights reserved.
*
* THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
* THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
* ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
* RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
* ACCORDING TO RECIPIENTS CHOICE.
*
* The OpenModelica software and the Open Source Modelica
* Consortium (OSMC) Public License (OSMC-PL) are obtained
* from OSMC, either from the above address,
* from the URLs: http:www.ida.liu.se/projects/OpenModelica or
* http:www.openmodelica.org, and in the OpenModelica distribution.
* GNU version 3 is obtained from: http:www.gnu.org/copyleft/gpl.html.
*
* This program is distributed WITHOUT ANY WARRANTY; without
* even the implied warranty of  MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
* IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
*
* See the full OSMC Public License conditions for more details.
*
*/ =#

module BackendDAE

using MetaModelica
using ExportAll

#= Predeclaration of mutable recursive types =#
@UniontypeDecl BackendDAE
@UniontypeDecl EqSystem
@UniontypeDecl SubClock
@UniontypeDecl BaseClockPartitionKind
@UniontypeDecl Shared
@UniontypeDecl InlineData
@UniontypeDecl BasePartition
@UniontypeDecl SubPartition
@UniontypeDecl PartitionsInfo
@UniontypeDecl ExtraInfo
@UniontypeDecl BackendDAEType
@UniontypeDecl DataReconciliationData
@UniontypeDecl Variables
@UniontypeDecl CrefIndex
@UniontypeDecl VariableArray
@UniontypeDecl Var
@UniontypeDecl VarKind
@UniontypeDecl TearingSelect
@UniontypeDecl EquationKind
@UniontypeDecl EvaluationStages
@UniontypeDecl EquationAttributes
@UniontypeDecl Equation
@UniontypeDecl WhenEquation
@UniontypeDecl WhenOperator
@UniontypeDecl ExternalObjectClass
@UniontypeDecl Matching
@UniontypeDecl IndexReduction
@UniontypeDecl EquationConstraints
@UniontypeDecl StateOrder
@UniontypeDecl StrongComponent
@UniontypeDecl TearingSet
@UniontypeDecl InnerEquation
@UniontypeDecl StateSet
@UniontypeDecl EventInfo
@UniontypeDecl ZeroCrossingSet
@UniontypeDecl ZeroCrossing
@UniontypeDecl TimeEvent
@UniontypeDecl Solvability
@UniontypeDecl IndexType
@UniontypeDecl JacobianType
@UniontypeDecl Jacobian
@UniontypeDecl DifferentiateInputData
@UniontypeDecl DifferentiationType
@UniontypeDecl CompInfo
@UniontypeDecl BackendDAEModeData

import Absyn
import DAE
import DoubleEnded
import SCode
import ZeroCrossings

const EquationArray = Array #= Let's use the Julia array instead ExpandableArray =#

Type = .DAE.Type  #= Once we are in BackendDAE, the Type can be only basic types or enumeration.
We cannot do this in DAE because functions may contain many more types.
adrpo: yes we can, we just simplify the DAE.Type, see Types.simplifyType =#

#= THE LOWERED DAE consist of variables and equations. The variables are split into
two lists, one for unknown variables states and algebraic and one for known variables
constants and parameters.
The equations are also split into two lists, one with simple equations, a=b, a-b=0, etc., that
are removed from  the set of equations to speed up calculations. =#
@Uniontype BackendDAE begin
  @Record DAE begin

    eqs::EqSystems
    shared::Shared
  end
end

EqSystems = List

#= An independent system of equations (and their corresponding variables) =#
@Uniontype EqSystem begin
  @Record EQSYSTEM begin

    orderedVars #= ordered Variables, only states and alg. vars =#::Variables
    orderedEqs #= ordered Equations =#::EquationArray
    m::Option{AdjacencyMatrix}
    mT::Option{AdjacencyMatrixT}
    mapping #= current type of adjacency matrix, boolean is true if scalar =#::Option{AdjacencyMatrixMapping}
    matching::Matching
    stateSets #= the state sets of the system =#::StateSets
    partitionKind::BaseClockPartitionKind
    removedEqs #= these are equations that cannot solve for a variable.
    e.g. assertions, external function calls, algorithm sections without effect =#::EquationArray
  end
end

#=Clock removed for now=#
@Uniontype SubClock begin
  @Record SUBCLOCK begin
  end

  @Record INFERED_SUBCLOCK begin
  end
end

const DEFAULT_SUBCLOCK = "TODO. NOT SUPPORTED"

@Uniontype BaseClockPartitionKind begin
  @Record UNKNOWN_PARTITION begin

  end

  @Record CLOCKED_PARTITION begin

    subPartIdx::Integer
  end

  @Record CONTINUOUS_TIME_PARTITION begin

  end

  @Record UNSPECIFIED_PARTITION begin

  end
end

#= Data shared for all equation-systems =#
@Uniontype Shared begin
  @Record SHARED begin

    globalKnownVars #= variables only depending on parameters and constants [TODO: move stuff (like inputs) to localKnownVars] =#::Variables
    localKnownVars #= variables only depending on locally constant variables in the simulation step, i.e. states, input variables =#::Variables
    externalObjects #= External object variables =#::Variables
    aliasVars #= Data originating from removed simple equations needed to build
    variables' lookup table (in C output).
    In that way, double buffering of variables in pre()-buffer, extrapolation
    buffer and results caching, etc., is avoided, but in C-code output all the
    data about variables' names, comments, units, etc. is preserved as well as
    pointer to their values (trajectories). =#::Variables
    initialEqs #= Initial equations =#::EquationArray
    removedEqs #= these are equations that cannot solve for a variable. for example assertions, external function calls, algorithm sections without effect =#::EquationArray
    constraints #= constraints (Optimica extension) =#::List{.DAE.Constraint}
    classAttrs #= class attributes (Optimica extension) =#::List{.DAE.ClassAttributes}
    cache #=TODO: Use the Julia cache here=#
    graph #=TODO use something else then the FCore crap=#
    functionTree #= functions for Backend =# #=Let's skip this part people=#
    eventInfo #= eventInfo =#::EventInfo
    extObjClasses #= classes of external objects, contains constructor & destructor =#::ExternalObjectClasses
    backendDAEType #= indicate for what the BackendDAE is used =#::BackendDAEType
    symjacs #= Symbolic Jacobians =#::SymbolicJacobians
    info #= contains extra info that we send around like the model name =#::ExtraInfo
    partitionsInfo::PartitionsInfo
    daeModeData #= DAEMode Data =#::BackendDAEModeData
    dataReconciliationData::Option{DataReconciliationData}
  end
end

@Uniontype InlineData begin
  @Record INLINE_DATA begin

    inlineSystems::EqSystems
    knownVariables::Variables
  end
end

@Uniontype BasePartition begin
  @Record BASE_PARTITION begin

    clock::.DAE.ClockKind
    nSubClocks::Integer
  end
end

@Uniontype SubPartition begin
  @Record SUB_PARTITION begin

    clock::SubClock
    holdEvents::Bool
    prevVars::List{.DAE.ComponentRef}
  end
end

@Uniontype PartitionsInfo begin
  @Record PARTITIONS_INFO begin

    basePartitions::Array{BasePartition}
    subPartitions::Array{SubPartition}
  end
end

#= extra information that we should send around with the DAE =#
@Uniontype ExtraInfo begin
  @Record EXTRA_INFO begin

    description #= the model description string =#::String
    fileNamePrefix #= the model name to be used in the dumps =#::String
  end
end

#= BackendDAEType to indicate different types of BackendDAEs.
For example for simulation, initialization, Jacobian, algebraic loops etc. =#
@Uniontype BackendDAEType begin
  @Record SIMULATION begin

  end

  @Record JACOBIAN begin

  end

  @Record ALGEQSYSTEM begin

  end

  @Record ARRAYSYSTEM begin

  end

  @Record PARAMETERSYSTEM begin

  end

  @Record INITIALSYSTEM begin

  end

  @Record INLINESYSTEM begin

  end

  @Record DAEMODESYSTEM begin

  end
end

@Uniontype DataReconciliationData begin
  @Record DATA_RECON begin

    symbolicJacobian #= SET_S w.r.t ... =#::Jacobian
    setcVars #= setc solved vars =#::Variables
    datareconinputs::Variables
    #=  ... maybe more DATA for the code generation
    =#
  end
end

#=
=#
#=   variables and equations definition
=#
#=
=#

@Uniontype Variables begin
  @Record VARIABLES begin

    crefIndices #= HashTB, cref->indx =#::Array{List{CrefIndex}}
    varArr #= Array of variables =#::VariableArray
    bucketSize #= bucket size =#::Integer
    numberOfVars #= no. of vars =#::Integer
  end
end

#= Component Reference Index =#
@Uniontype CrefIndex begin
  @Record CREFINDEX begin

    cref::.DAE.ComponentRef
    index::Integer
  end
end

#= array of Equations are expandable, to amortize the cost of adding
equations in a more efficient manner =#
@Uniontype VariableArray begin
  @Record VARIABLE_ARRAY begin

    numberOfElements #= no. elements =#::Integer
    varOptArr::Array{Option{Var}}
  end
end


#= variables =#
@Uniontype Var begin
  @Record VAR begin

    varName #= variable name =#::.DAE.ComponentRef
    varKind #= kind of variable =#::VarKind
    varDirection #= input, output or bidirectional =#::.DAE.VarDirection
    varParallelism #= parallelism of the variable. parglobal, parlocal or non-parallel =#::.DAE.VarParallelism
    varType #= built-in type or enumeration =#::Type
    bindExp #= Binding expression e.g. for parameters =#::Option{.DAE.Exp}
    tplExp #= Variable is part of a tuple. Needed for the globalKnownVars and localKnownVars =#::Option{.DAE.Exp}
    arryDim #= array dimensions of non-expanded var =#::.DAE.InstDims
    source #= origin of variable =#::.DAE.ElementSource
    values #= values on built-in attributes =#::Option{.DAE.VariableAttributes}
    tearingSelectOption #= value for TearingSelect =#::Option{TearingSelect}
    hideResult #= expression from the hideResult annotation =#::.DAE.Exp
    comment #= this contains the comment and annotation from Absyn =#::Option{SCode.Comment}
    connectorType #= flow, stream, unspecified or not connector. =#::.DAE.ConnectorType
    innerOuter #= inner, outer, inner outer or unspecified =#::.DAE.VarInnerOuter
    unreplaceable #= indicates if it is allowed to replace this variable =#::Bool
  end
end

#= variable kind =#
@Uniontype VarKind begin
  @Record VARIABLE begin

  end

  @Record STATE begin

    index #= how often this states was differentiated =#::Integer
    derName #= the name of the derivative =#::Option{.DAE.ComponentRef}
    natural #= false if it was forced by StateSelect.always or StateSelect.prefer or generated by index reduction =#::Bool
  end

  @Record STATE_DER begin

  end

  @Record DUMMY_DER begin

  end

  @Record DUMMY_STATE begin

  end

  @Record CLOCKED_STATE begin

    previousName #= the name of the previous variable =#::.DAE.ComponentRef
    isStartFixed #= is fixed at first clock tick =#::Bool
  end

  @Record DISCRETE begin

  end

  @Record PARAM begin

  end

  @Record CONST begin

  end

  @Record EXTOBJ begin

    fullClassName::Absyn.Path
  end

  @Record JAC_VAR begin

  end

  @Record JAC_DIFF_VAR begin

  end

  @Record SEED_VAR begin

  end

  @Record OPT_CONSTR begin

  end

  @Record OPT_FCONSTR begin

  end

  @Record OPT_INPUT_WITH_DER begin

  end

  @Record OPT_INPUT_DER begin

  end

  @Record OPT_TGRID begin

  end

  @Record OPT_LOOP_INPUT begin

    replaceExp::.DAE.ComponentRef
  end

  @Record ALG_STATE begin

  end

  @Record ALG_STATE_OLD begin

  end

  @Record DAE_RESIDUAL_VAR begin

  end

  @Record DAE_AUX_VAR begin

  end

  @Record LOOP_ITERATION begin

  end

  @Record LOOP_SOLVED begin

  end
end

@Uniontype TearingSelect begin
  @Record NEVER begin

  end

  @Record AVOID begin

  end

  @Record DEFAULT begin

  end

  @Record PREFER begin

  end

  @Record ALWAYS begin

  end
end

const WHENCLK_PRREFIX = "whenclk"::String

#= equation kind =#
@Uniontype EquationKind begin
  @Record BINDING_EQUATION begin

  end

  @Record DYNAMIC_EQUATION begin

  end

  @Record INITIAL_EQUATION begin

  end

  @Record CLOCKED_EQUATION begin

    clk::Integer
  end

  @Record DISCRETE_EQUATION begin

  end

  @Record AUX_EQUATION begin

  end

  @Record UNKNOWN_EQUATION_KIND begin

  end
end

#= evaluation stages =#
@Uniontype EvaluationStages begin
  @Record EVALUATION_STAGES begin

    dynamicEval::Bool
    algebraicEval::Bool
    zerocrossEval::Bool
    discreteEval::Bool
  end
end

const defaultEvalStages = EVALUATION_STAGES(false, false, false, false)::EvaluationStages

@Uniontype EquationAttributes begin
  @Record EQUATION_ATTRIBUTES begin

    differentiated #= true if the equation was differentiated, and should not be differentiated again to avoid equal equations =#::Bool
    kind::EquationKind
    evalStages::EvaluationStages
  end
end

const EQ_ATTR_DEFAULT_DYNAMIC = EQUATION_ATTRIBUTES(false, DYNAMIC_EQUATION(), defaultEvalStages)::EquationAttributes

const EQ_ATTR_DEFAULT_BINDING = EQUATION_ATTRIBUTES(false, BINDING_EQUATION(), defaultEvalStages)::EquationAttributes

const EQ_ATTR_DEFAULT_INITIAL = EQUATION_ATTRIBUTES(false, INITIAL_EQUATION(), defaultEvalStages)::EquationAttributes

const EQ_ATTR_DEFAULT_DISCRETE = EQUATION_ATTRIBUTES(false, DISCRETE_EQUATION(), defaultEvalStages)::EquationAttributes

const EQ_ATTR_DEFAULT_AUX = EQUATION_ATTRIBUTES(false, AUX_EQUATION(), defaultEvalStages)::EquationAttributes

const EQ_ATTR_DEFAULT_UNKNOWN = EQUATION_ATTRIBUTES(false, UNKNOWN_EQUATION_KIND(), defaultEvalStages)::EquationAttributes

@Uniontype Equation begin
  @Record EQUATION begin

    exp::.DAE.Exp
    scalar::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record ARRAY_EQUATION begin

    dimSize #= dimension sizes =#::List{Integer}
    left #= lhs =#::.DAE.Exp
    right #= rhs =#::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
    recordSize #= NONE() if not a record =#::Option{Integer}
  end

  @Record SOLVED_EQUATION begin

    componentRef::.DAE.ComponentRef
    exp::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record RESIDUAL_EQUATION begin

    exp #= not present from FrontEnd =#::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record ALGORITHM begin

    size #= size of equation =#::Integer
    alg::.DAE.Algorithm
    source #= origin of algorithm =#::.DAE.ElementSource
    expand #= this algorithm was translated from an equation. we should not expand array crefs! =#::.DAE.Expand
    attr::EquationAttributes
  end

  @Record WHEN_EQUATION begin

    size #= size of equation =#::Integer
    whenEquation::WhenEquation
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record COMPLEX_EQUATION begin

    size #= size of equation =#::Integer
    left #= lhs =#::.DAE.Exp
    right #= rhs =#::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record IF_EQUATION begin

    conditions #= Condition =#::List{.DAE.Exp}
    eqnstrue #= Equations of true branch =#::List{List{Equation}}
    eqnsfalse #= Equations of false branch =#::List{Equation}
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record FOR_EQUATION begin

    iter #= the iterator variable =#::.DAE.Exp
    start #= start of iteration =#::.DAE.Exp
    stop #= end of iteration =#::.DAE.Exp
    body #= iterated equation =#::Equation
    source #= origin of equation =#::.DAE.ElementSource
    attr::EquationAttributes
  end

  @Record DUMMY_EQUATION begin

  end
end

@Uniontype WhenEquation begin
  @Record WHEN_STMTS begin

    condition #= the when-condition =#::.DAE.Exp
    whenStmtLst::List{WhenOperator}
    elsewhenPart #= elsewhen equation with the same cref on the left hand side. =#::Option{WhenEquation}
  end
end

@Uniontype WhenOperator begin
  @Record ASSIGN begin

    left #= left hand side of equation =#::.DAE.Exp
    right #= right hand side of equation =#::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
  end

  @Record REINIT begin

    stateVar #= State variable to reinit =#::.DAE.ComponentRef
    value #= Value after reinit =#::.DAE.Exp
    source #= origin of equation =#::.DAE.ElementSource
  end

  @Record ASSERT begin

    condition::.DAE.Exp
    message::.DAE.Exp
    level::.DAE.Exp
    source #= the origin of the component/equation/algorithm =#::.DAE.ElementSource
  end

  @Record TERMINATE begin

    message::.DAE.Exp
    source #= the origin of the component/equation/algorithm =#::.DAE.ElementSource
  end

  @Record NORETCALL begin

    exp::.DAE.Exp
    source #= the origin of the component/equation/algorithm =#::.DAE.ElementSource
  end
end

ExternalObjectClasses = List  #= classes of external objects stored in list =#

#= class of external objects =#
@Uniontype ExternalObjectClass begin
  @Record EXTOBJCLASS begin

    path #= className of external object =#::Absyn.Path
    source #= origin of equation =#::.DAE.ElementSource
  end
end

#=
=#
#=   Matching, strong components and StateSets
=#
#=
=#

@Uniontype Matching begin
  @Record NO_MATCHING begin

  end

  @Record MATCHING begin

    ass1 #= ass[varindx]=eqnindx =#::Array{Integer}
    ass2 #= ass[eqnindx]=varindx =#::Array{Integer}
    comps::StrongComponents
  end
end

@Uniontype IndexReduction begin
  @Record INDEX_REDUCTION begin

  end

  @Record NO_INDEX_REDUCTION begin

  end
end

@Uniontype EquationConstraints begin
  @Record ALLOW_UNDERCONSTRAINED begin

  end

  @Record EXACT begin

  end
end

MatchingOptions = Tuple

StructurallySingularSystemHandlerArg = Tuple  #= StateOrder,ConstraintEqns,Eqn->EqnsIndxes,EqnIndex->Eqns,NrOfEqnsbeforeIndexReduction =#

ConstraintEquations = Array

@Uniontype StateOrder begin
  @Record STATEORDER begin

    hashTable #= x -> dx =#
    invHashTable #= dx -> {x,y,z} =#
  end

  @Record NOSTATEORDER begin

  end
end

StrongComponents = List  #= Order of the equations the have to be solved =#

@Uniontype StrongComponent begin
  @Record SINGLEEQUATION begin

    eqn::Integer
    var::Integer
  end

  @Record EQUATIONSYSTEM begin

    eqns::List{Integer}
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
    jac::Jacobian
    jacType::JacobianType
    mixedSystem #= true for system that discrete dependencies to the iteration variables =#::Bool
  end

  @Record SINGLEARRAY begin

    eqn::Integer
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
  end

  @Record SINGLEALGORITHM begin

    eqn::Integer
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
  end

  @Record SINGLECOMPLEXEQUATION begin

    eqn::Integer
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
  end

  @Record SINGLEWHENEQUATION begin

    eqn::Integer
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
  end

  @Record SINGLEIFEQUATION begin

    eqn::Integer
    vars #= be careful with states, this are solved for der(x) =#::List{Integer}
  end

  @Record TORNSYSTEM begin

    strictTearingSet::TearingSet
    casualTearingSet::Option{TearingSet}
    linear::Bool
    mixedSystem #= true for system that discrete dependencies to the iteration variables =#::Bool
  end
end

@Uniontype TearingSet begin
  @Record TEARINGSET begin

    tearingvars::List{Integer}
    residualequations::List{Integer}
    innerEquations #= list of matched equations and variables; these will be solved explicitly in the given order =#::InnerEquations
    jac::Jacobian
  end
end

InnerEquations = List

@Uniontype InnerEquation begin
  @Record INNEREQUATION begin

    eqn::Integer
    vars::List{Integer}
  end

  @Record INNEREQUATIONCONSTRAINTS begin

    eqn::Integer
    vars::List{Integer}
    cons::Constraints
  end
end

StateSets = List  #= List of StateSets =#

@Uniontype StateSet begin
  @Record STATESET begin

    index::Integer
    rang::Integer
    #=  how many states are needed?
    =#
    state::List{.DAE.ComponentRef}
    crA #= set.x=A*states =#::.DAE.ComponentRef
    varA::List{Var}
    #= the jacobian matrix entries
    =#
    statescandidates::List{Var}
    #= all state candidates
    =#
    ovars::List{Var}
    #= other variables to solve the eqns
    =#
    eqns::List{Equation}
    #= the constraint equations
    =#
    oeqns::List{Equation}
    #= other equations to solve the eqns
    =#
    crJ::.DAE.ComponentRef
    #=  the jac vector
    =#
    varJ::List{Var}
    jacobian::Jacobian
  end
end

#=
=#
#=  event info and stuff
=#
#=
=#

@Uniontype EventInfo begin
  @Record EVENT_INFO begin

    timeEvents #= stores all information related to time events =#::List{TimeEvent}
    zeroCrossings #= list of zero crossing conditions =#::ZeroCrossingSet
    relations #= list of zero crossing function as before =#::DoubleEnded.MutableList{ZeroCrossing}
    samples #= [deprecated] list of sample as before, only used by cpp runtime (TODO: REMOVE ME) =#::ZeroCrossingSet
    numberMathEvents #= stores the number of math function that trigger events e.g. floor, ceil, integer, ... =#::Integer
  end
end

@Uniontype ZeroCrossingSet begin
  @Record ZERO_CROSSING_SET begin

    zc::DoubleEnded.MutableList{ZeroCrossing}
    tree::Array{ZeroCrossings.Tree}
  end
end

@Uniontype ZeroCrossing begin
  @Record ZERO_CROSSING begin

    relation_ #= function =#::.DAE.Exp
    occurEquLst #= list of equations where the function occurs =#::List{Integer}
  end
end

@Uniontype TimeEvent begin
  @Record SIMPLE_TIME_EVENT begin

  end

  @Record SAMPLE_TIME_EVENT begin

    index #= unique sample index =#::Integer
    startExp::.DAE.Exp
    intervalExp::.DAE.Exp
  end
end

#=
=#
#=  AdjacencyMatrixes
=#
#=
=#

AdjacencyMatrixElementEntry = Integer
AdjacencyMatrixElement = List
AdjacencyMatrix = Array  #= array<list<Integer>> =#
AdjacencyMatrixT = AdjacencyMatrix  #= a list of equation indices (1..n), one for each variable. Equations that -only-
contain the state variable and not the derivative have a negative index. =#
AdjacencyMatrixMapping = Tuple  #= a mapping for adjacency matrices that contains:
array<list<Integer>>: array index -> scalar index list
array<Integer>      : scalar index -> array index (not unique)
IndexType           : the occurence condition type for the current adjacency matrix
Boolean             : true if scalar
Boolean             : true if analytical to structural singularity processing has already been done =#

AdjacencyMatrixElementEnhancedEntry = Tuple
AdjacencyMatrixElementEnhanced = List
AdjacencyMatrixEnhanced = Array
AdjacencyMatrixTEnhanced = AdjacencyMatrixEnhanced

@Uniontype Solvability begin
  @Record SOLVABILITY_SOLVED begin

  end

  @Record SOLVABILITY_CONSTONE begin

  end

  @Record SOLVABILITY_CONST begin

    b #= false if the constant is almost zero (<1e-6) =#::Bool
  end

  @Record SOLVABILITY_PARAMETER begin

    b #= false if the partial derivative is zero =#::Bool
  end

  @Record SOLVABILITY_LINEAR begin

    b #= false if the partial derivative is zero =#::Bool
  end

  @Record SOLVABILITY_NONLINEAR begin

  end

  @Record SOLVABILITY_UNSOLVABLE begin

  end

  @Record SOLVABILITY_SOLVABLE begin

  end
end

Constraints = List  #= Constraints on the solvability of the (casual) tearing set; needed for proper Dynamic Tearing =#

@Uniontype IndexType begin
  @Record ABSOLUTE begin

  end

  @Record NORMAL begin

  end

  @Record SOLVABLE begin

  end

  @Record BASECLOCK_IDX begin

  end

  @Record SUBCLOCK_IDX begin

  end

  @Record SPARSE begin

  end
end

#=
=#
#=  Jacobian stuff
=#
#=
=#

@Uniontype JacobianType begin
  @Record JAC_CONSTANT begin

  end

  @Record JAC_LINEAR begin

  end

  @Record JAC_NONLINEAR begin

  end

  @Record JAC_GENERIC begin

  end

  @Record JAC_NO_ANALYTIC begin

  end
end

const SymbolicJacobianAIndex = 1::Integer

const SymbolicJacobianBIndex = 2::Integer

const SymbolicJacobianCIndex = 3::Integer

const SymbolicJacobianDIndex = 4::Integer

const derivativeNamePrefix = "DERAlias"::String

const partialDerivativeNamePrefix = "pDER"::String

const functionDerivativeNamePrefix = "funDER"::String

const outputAliasPrefix = "outputAlias_"::String

const optimizationMayerTermName = "OMCobjectMayerTerm"::String

const optimizationLagrangeTermName = "OMCobjectLagrangeTerm"::String

const symSolverDT = "__OMC_DT"::String

const homotopyLambda = "__HOM_LAMBDA"::String
FullJacobian = Option

@Uniontype Jacobian begin
  @Record FULL_JACOBIAN begin

    jacobian::FullJacobian
  end

  @Record GENERIC_JACOBIAN begin

    jacobian::Option{SymbolicJacobian}
    sparsePattern::SparsePattern
    coloring::SparseColoring
  end

  @Record EMPTY_JACOBIAN begin

  end
end

SymbolicJacobians = List

SymbolicJacobian = Tuple
#=  symbolic equation system
=#
#=  Matrix name
=#
#=  diff vars (independent vars)
=#
#=  diffed vars (residual vars)
=#
#=  all diffed vars (residual vars + dependent vars)
=#
#=  original dependent variables
=#
SparsePatternCref = Tuple
SparsePatternCrefs = List

SparsePattern = Tuple
#=  column-wise sparse pattern
=#
#=  row-wise sparse pattern
=#
#=  diff vars (independent vars) of associated jacobian
=#
#=  diffed vars (residual vars) of associated jacobian
=#
#=  nonZeroElements
=#

const emptySparsePattern = (nil, nil, (nil, nil), 0)::SparsePattern

SparseColoring = List
#=  colouring
=#
#= /*
Type only for transformation from analytical to structural singularity.
*/ =#
LinearIntegerJacobianRow = List
#=  Actual jacobian entries sparse, <column, value>
=#
LinearIntegerJacobianRhs = Array
#=  RHS-Exp for full pivot algorithm. Replacement for eliminated equation.
=#
LinearIntegerJacobianIndices = Array
#=  Index tuple <array, scalar> for equations
=#
#= /*
The full linear integer matrix
- additional boolean array to track which rows have been changed
- additional boolean array to track which variables are matched to the equations
*/ =#
LinearIntegerJacobian = Tuple

@Uniontype DifferentiateInputData begin
  @Record DIFFINPUTDATA begin

    independenentVars::Option{Variables}
    #=  Independent variables
    =#
    dependenentVars::Option{Variables}
    #=  Dependent variables
    =#
    knownVars::Option{Variables}
    #=  known variables (e.g. parameter, constants, ...)
    =#
    allVars::Option{Variables}
    #=  all variables
    =#
    controlVars::List{Var}
    #=  variables to save control vars of for algorithm
    =#
    diffCrefs::List{.DAE.ComponentRef}
    #=  all crefs to differentiate, needed for generic gradient
    =#
    matrixName::Option{String}
    #=  name to create temporary vars, needed for generic gradient
    =#
    diffedFunctions #=TODO: Use something sane here=#
    #=  current functions, to prevent recursive differentiation
    =#
  end
end

const emptyInputData = DIFFINPUTDATA(NONE(), NONE(), NONE(), NONE(), nil, nil, NONE(), nothing #=TODO: Fixme=#)::DifferentiateInputData

DifferentiateInputArguments = Tuple

#= Define the behaviour of differentiation method for (e.g. index reduction, ...) =#
@Uniontype DifferentiationType begin
  @Record DIFFERENTIATION_TIME begin

  end

  @Record SIMPLE_DIFFERENTIATION begin

  end

  @Record DIFFERENTIATION_FUNCTION begin

  end

  @Record DIFF_FULL_JACOBIAN begin

  end

  @Record GENERIC_GRADIENT begin

  end
end

#= types to count operations for the components =#
@Uniontype CompInfo begin
  @Record COUNTER begin

    #=  single equation
    =#
    comp::StrongComponent
    numAdds::Integer
    numMul::Integer
    numDiv::Integer
    numTrig::Integer
    numRelations::Integer
    numLog::Integer
    #=  logical operations
    =#
    numOth::Integer
    #=  pow,...
    =#
    funcCalls::Integer
  end

  @Record SYSTEM begin

    #= linear system of equations
    =#
    comp::StrongComponent
    allOperations::CompInfo
    size::Integer
    density::Float
  end

  @Record TORN_ANALYSE begin

    #= torn system of equations
    =#
    comp::StrongComponent
    tornEqs::CompInfo
    otherEqs::CompInfo
    tornSize::Integer
  end

  @Record NO_COMP begin

    #=  assert...
    =#
    numAdds::Integer
    numMul::Integer
    numDiv::Integer
    numTrig::Integer
    numRelations::Integer
    numLog::Integer
    #=  logical operations
    =#
    numOth::Integer
    #=  pow,...
    =#
    funcCalls::Integer
  end
end

@Uniontype BackendDAEModeData begin
  @Record BDAE_MODE_DATA begin

    stateVars::List{Var}
    algStateVars::List{Var}
    numResVars::Integer
    modelVars::Option{Variables}
  end
end

const emptyDAEModeData = BDAE_MODE_DATA(nil, nil, 0, NONE())::BackendDAEModeData

#= So that we can use wildcard imports and named imports when they do occur. Not good Julia practice =#
@exportAll()

end
