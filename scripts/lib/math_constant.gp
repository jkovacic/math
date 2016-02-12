/*
 * A PARI/GP script that reproduces mathematical constants
 * that are defined in 'lib/util/math_constant.h', in
 * arbitrary precision. PARI/GP default is 38 characters.
 *
 * From a shell, run the script as:
 *   gp </path/to/math_constant.gp
 */


mathConst() =
{
  print("Pi = ", Pi );
  print("1/Pi = ", 1/Pi );
  print("sqrt(Pi) = ", sqrt(Pi) );
  print("sqrt(1/Pi) = ", sqrt(1/Pi) );
  print("sqrt(2/Pi) = ", sqrt(2/Pi) );
  print("sqrt(1/(2*Pi)) = ", sqrt(1/(2*Pi)) );
  print("log(sqrt(2*Pi)) = ", log(sqrt(2*Pi)) );

  print("sqrt(2) = ", sqrt(2) );
  print("sqrt(2)/2 = ", sqrt(2)/2 );
}


mathConst()
