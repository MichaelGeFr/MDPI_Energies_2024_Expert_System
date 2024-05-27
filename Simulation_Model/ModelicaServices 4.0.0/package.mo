package ModelicaServices "ModelicaServices (Dymola implementation) - Models and functions used in the Modelica Standard Library requiring a tool specific implementation"
extends Modelica.Icons.Package;

package Machine

  final constant Real eps=1.e-15 "Biggest number such that 1.0 + eps = 1.0";

  final constant Real small=1.e-60
    "Smallest number such that small and -small are representable on the machine";

  final constant Real inf=1.e+60
    "Biggest Real number such that inf and -inf are representable on the machine";
  annotation (Documentation(info="<html>
<p>
Package in which processor specific constants are defined that are needed
by numerical algorithms. Typically these constants are not directly used,
but indirectly via the alias definition in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.
</p>
</html>"));
end Machine;
annotation (
  Protection(access=Access.hide),
  preferredView="info",
  version="4.0.0",
  versionBuild=1,
  versionDate="2020-01-15",
  dateModified = "2023-01-25 10:08:00Z",
  uses(Modelica(version="4.0.0")),
  conversion(
    noneFromVersion="1.0",
    noneFromVersion="1.1",
    noneFromVersion="1.2",
    noneFromVersion="3.2.1",
 noneFromVersion="3.2.3"),
  Documentation(info="<html>
<p>
This package contains a set of functions and models to be used in the
Modelica Standard Library that requires a tool specific implementation.
These are:
</p>

<ul>
<li> <a href=\"modelica://ModelicaServices.Animation.Shape\">Animation.Shape</a>
     provides a 3-dim. visualization of elementary
     mechanical objects. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Animation.Surface\">Animation.Surface</a>
     provides a 3-dim. visualization of
     moveable parameterized surface. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Animation.Vector\">Animation.Vector</a>
     provides a 3-dim. visualization of a vector objects. It is used in
<a href=\"modelica://Modelica.Mechanics.MultiBody.Visualizers.Advanced.Vector\">Modelica.Mechanics.MultiBody.Visualizers.Advanced.Vector</a>
     via inheritance.</li>
	 
<li> <a href=\"modelica://ModelicaServices.ExternalReferences.loadResource\">ExternalReferences.loadResource</a>
     provides a function to return the absolute path name of an URI or a local file name. It is used in
<a href=\"modelica://Modelica.Utilities.Files.loadResource\">Modelica.Utilities.Files.loadResource</a>
     via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Machine\">Machine</a>
     provides a package of machine constants. It is used in
<a href=\"modelica://Modelica.Constants\">Modelica.Constants</a>.</li>

<li> <a href=\"modelica://ModelicaServices.System.exit\">System.exit</a> provides a function to terminate the execution of the Modelica environment. It is used in <a href=\"modelica://Modelica.Utilities.System.exit\">Modelica.Utilities.System.exit</a> via inheritance.</li>

<li> <a href=\"modelica://ModelicaServices.Types.SolverMethod\">Types.SolverMethod</a>
     provides a string defining the integration method to solve differential equations in
     a clocked discretized continuous-time partition (see Modelica 3.3 language specification).
     It is not yet used in the Modelica Standard Library, but in the Modelica_Synchronous library
     that provides convenience blocks for the clock operators of Modelica version &ge; 3.3.</li>
</ul>

<p>
This is the Dymola implementation.
</p>

<p>
Original version
<strong>Licensed by the Modelica Association under the 3-Clause BSD License</strong><br>
Copyright &copy; 2009-2020, Modelica Association and contributors.
</p>
<p>
<strong>Modifications licensed by Dassault Syst&egrave;mes, copyright &copy; 2009-2020.</strong>
</p>

</html>"));
end ModelicaServices;
