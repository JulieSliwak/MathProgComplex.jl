param NONE 			symbolic := 'NONE';
param VAR_TYPE 	symbolic := 'VAR_TYPE';
param MONO_DEF 	symbolic := 'MONO_DEF';
param MONO		 	symbolic := 'MONO';
param QUAD 			symbolic := 'QUAD';
param LIN 			symbolic := 'LIN';
param CONST			symbolic := 'CONST';
param OBJ 			symbolic := 'OBJ';
param LB 				symbolic := 'LB';
param UB 				symbolic := 'UB';

/* Variable types */
param REAL symbolic := 'REAL';
param BOOL symbolic := 'BOOL';

/* Preconditioning types */
param SQRT 				symbolic := 'SQRT';


/* Problem description */
set KEYS dimen 4;
param RIGHT{KEYS} symbolic;
param LEFT {KEYS} symbolic;

/* Preconditioning of some constraints */
set PRECOND_CONSTRAINTS dimen 2;



set REAL_VARIABLES := setof{(VAR_TYPE, REAL, name, NONE) in KEYS}name;
set BINARY_VARIABLES := setof{(VAR_TYPE, BOOL, name, NONE) in KEYS}name;

var x{REAL_VARIABLES};
var y{BINARY_VARIABLES} binary;



set PRECOND_CSTRS := setof{(cstrname, any) in PRECOND_CONSTRAINTS}cstrname;
set SQRT_PRECOND_CSTR := setof{(cstrname, SQRT) in PRECOND_CONSTRAINTS}cstrname;

set CONSTRAINTS :=
	setof{(LB, name, NONE, NONE) in KEYS}name
	union
	setof{(UB, name, NONE, NONE) in KEYS}name;

set CONSTRAINTS_LB := setof{(LB, name, NONE, NONE) in KEYS}name;
set CONSTRAINTS_UB := setof{(UB, name, NONE, NONE) in KEYS}name;



minimize CRITERION:
  +sum{(MONO, OBJ, monname, NONE) in KEYS}(
		+LEFT [MONO, OBJ, monname, NONE] * prod{(MONO_DEF, monname, var0, NONE) in KEYS}(
			(if var0 in REAL_VARIABLES then x[var0] else y[var0])^LEFT [MONO_DEF, monname, var0, NONE]
		)
	)
	+sum{(QUAD, OBJ, var1, var2) in KEYS}(
		+LEFT [QUAD, OBJ, var1, var2] * (if var1 in REAL_VARIABLES then x[var1] else y[var1]) * (if var2 in REAL_VARIABLES then x[var2] else y[var2])
	)
	+sum{(LIN, OBJ, var1, NONE) in KEYS}(
		+LEFT [LIN, OBJ, var1, NONE] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
	)
	+sum{(LIN, OBJ, NONE, var1) in KEYS}(
		+LEFT [LIN, OBJ, NONE, var1] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
	)
 +(if (CONST, OBJ, NONE, NONE) in KEYS then LEFT[CONST, OBJ, NONE, NONE] else 0)
	;



subject to constraint{ctr in CONSTRAINTS diff SQRT_PRECOND_CSTR}:
	+(if ctr in CONSTRAINTS_LB then LEFT[LB, ctr, NONE, NONE] else -Infinity)
	<=
	+sum{(MONO, ctr, monname, NONE) in KEYS}(
		+LEFT [MONO, ctr, monname, NONE] * prod{(MONO_DEF, monname, var0, NONE) in KEYS}(
			(if var0 in REAL_VARIABLES then x[var0] else y[var0])^LEFT [MONO_DEF, monname, var0, NONE]
		)
	)
	+sum{(QUAD, ctr, var1, var2) in KEYS}(
		+LEFT [QUAD, ctr, var1, var2] * (if var1 in REAL_VARIABLES then x[var1] else y[var1]) * (if var2 in REAL_VARIABLES then x[var2] else y[var2])
	)
	+sum{(LIN, ctr, var1, NONE) in KEYS}(
		+LEFT [LIN, ctr, var1, NONE] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
	)
	+sum{(LIN, ctr, NONE, var1) in KEYS}(
		+LEFT [LIN, ctr, NONE, var1] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
	)
 	+(if (CONST, ctr, NONE, NONE) in KEYS then LEFT[CONST, ctr, NONE, NONE] else 0)
	<=
	+(if ctr in CONSTRAINTS_UB then LEFT[UB, ctr, NONE, NONE] else +Infinity)
	;


/* SQRT precond constraints */
subject to sqrt_constraint{ctr in SQRT_PRECOND_CSTR}:
	#+(if ctr in CONSTRAINTS_LB then LEFT[LB, ctr, NONE, NONE] else -Infinity)
	+(if ctr in CONSTRAINTS_LB then sqrt(LEFT[LB, ctr, NONE, NONE]) else -Infinity)
	<=
	sqrt(
		+1e-6
		+sum{(MONO, ctr, monname, NONE) in KEYS}(
			+LEFT [MONO, ctr, monname, NONE] * prod{(MONO_DEF, monname, var0, NONE) in KEYS}(
				(if var0 in REAL_VARIABLES then x[var0] else y[var0])^LEFT [MONO_DEF, monname, var0, NONE]
			)
		)
		+sum{(QUAD, ctr, var1, var2) in KEYS}(
			+LEFT [QUAD, ctr, var1, var2] * (if var1 in REAL_VARIABLES then x[var1] else y[var1]) * (if var2 in REAL_VARIABLES then x[var2] else y[var2])
		)
		+sum{(LIN, ctr, var1, NONE) in KEYS}(
			+LEFT [LIN, ctr, var1, NONE] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
		)
		+sum{(LIN, ctr, NONE, var1) in KEYS}(
			+LEFT [LIN, ctr, NONE, var1] * (if var1 in REAL_VARIABLES then x[var1] else y[var1])
		)
	 	+(if (CONST, ctr, NONE, NONE) in KEYS then LEFT[CONST, ctr, NONE, NONE] else 0)
	)
	<=
	+(if ctr in CONSTRAINTS_UB then sqrt(LEFT[UB, ctr, NONE, NONE]) else +Infinity)
	#+(if ctr in CONSTRAINTS_UB then LEFT[UB, ctr, NONE, NONE] else +Infinity)
	;
