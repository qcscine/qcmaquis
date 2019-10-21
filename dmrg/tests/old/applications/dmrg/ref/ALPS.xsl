<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<!--
   Copyright (c) 2003-2010 Matthias Troyer (troyer@ethz.ch)
    
   Distributed under the Boost Software License, Version 1.0.
   (See accompanying file LICENSE_1_0.txt or copy at
   http://www.boost.org/LICENSE_1_0.txt)
  -->

  <xsl:template match="/TEST">
    <html>
      <head/>
      <body>
        <H1>ALPS Job Test</H1>
          <xsl:for-each select="INPUT">
            <table border="0"> <tbody>
            <tr> <td> Test Input File </td> <td> <A href="{@file}"><xsl:value-of select="@file"/> </A> </td> </tr>
            </tbody> </table>
          </xsl:for-each>

          <xsl:for-each select="EXECUTED">
          <H2>Execution information</H2>
            <table border="1"><thead><tr><td><b>Executed from</b></td><td><b>to</b></td>
              <td><b>on</b></td></tr></thead><tbody>
              <tr>
                <td><xsl:for-each select="FROM"><xsl:apply-templates/></xsl:for-each></td>
                <td><xsl:for-each select="TO"><xsl:apply-templates/></xsl:for-each></td>
                <td><xsl:for-each select="MACHINE">
                      <xsl:for-each select="NAME"><xsl:apply-templates/></xsl:for-each>
                    </xsl:for-each></td>
              </tr>
            </tbody></table>
          </xsl:for-each>

          <xsl:for-each select="SCRIPT">
            <H2>Test Information</H2>
            <p>
            <H3>Executable</H3>
              <A href="{@file}"> <xsl:value-of select="@file"/> </A>
            </p>
          </xsl:for-each>

          <xsl:for-each select="SCRIPTINPUTFILES">
            <p>
            <table border="0"> <tbody>
            <xsl:for-each select="SCRIPTINPUTFILE">
                <tr> <td>Input File</td> <td> <A href="{@file}"> <xsl:value-of select="@file"/> </A> </td> </tr>
            </xsl:for-each>
            </tbody> </table>
            </p>
          </xsl:for-each>

          <xsl:for-each select="TESTINFOS">
          <H3>Test Settings</H3>
            <table border="1"> <tbody>
            <xsl:for-each select="TESTINFO">
              <tr> <td><xsl:value-of select="@name"/></td> <td> <xsl:apply-templates/> </td> </tr>
            </xsl:for-each>
            </tbody> </table>
          </xsl:for-each>

          <xsl:for-each select="REFERENCE">
            <H2>Reference Data</H2>
            <table border="0"> <tbody>
            <xsl:for-each select="FILE">
              <tr> <td> <A href="{@file}"> <xsl:value-of select="@file"/> </A> </td> </tr>
            </xsl:for-each>
            </tbody> </table>
          </xsl:for-each>

          <xsl:for-each select="OBSERVABLES">
            <H2>List Of Tested Observables</H2>
            <table border="0"> <tbody>
            <xsl:for-each select="OBSERVABLE">
              <tr> <td> <xsl:apply-templates/> </td> </tr>
            </xsl:for-each>
            </tbody> </table>
          </xsl:for-each>

          <H2>Test Results</H2>
          <xsl:for-each select="OUTPUT">
            <A href="{@file}"><xsl:value-of select="@file"/></A>
          </xsl:for-each>
          <xsl:for-each select="TESTTASK">
              <p>
              <table border="1">
                  <thead><tr><td><b>File</b></td><td><b>Result</b></td>
                      </tr></thead>
                  <tbody>
                    <tr>
                       <td><xsl:value-of select="@file"/></td>
                       <xsl:choose>
                         <xsl:when test="@status = 'PASS'">
                           <td bgcolor="#00ff00"><xsl:value-of select="@status"/></td>
                         </xsl:when>
                         <xsl:otherwise>
                           <td bgcolor="#ff0000"><xsl:value-of select="@status"/></td>
                         </xsl:otherwise>
                       </xsl:choose>
                    </tr>
                  </tbody>
              </table>
              </p>
              <p style='margin-left:40px;'>
              <table border="1">
                <xsl:for-each select="OBSERVABLES">
                  <thead><tr><td><b>Observable</b></td><td><b>|test-ref|</b></td><td><b>max. tolerance</b></td></tr>
                  </thead>
                  <tbody>
                    <xsl:for-each select="OBSERVABLE">
                      <tr> <td><xsl:value-of select="@name"/></td>
                           <td><xsl:value-of select="@diff"/></td>
                           <td><xsl:value-of select="@tol"/></td>
                      </tr>
                    </xsl:for-each>
                  </tbody>
                </xsl:for-each>
              </table>
              </p>
          </xsl:for-each>

          <xsl:for-each select="OARCHIVE">
              <H2> Compressed Output </H2>
              <A href="{@file}"> <xsl:value-of select="@file"/> </A>
          </xsl:for-each>
          
      </body>
    </html>
  </xsl:template>

  <xsl:template match="/JOB">
    <html>
      <head/>
      <body>
        <H1>ALPS Job Description File</H1>
        <xsl:for-each select="OUTPUT">
          <p>Output to new job file <A href="{@file}"></A></p>
        </xsl:for-each>
        <table border="1">
          <thead><tr><td><b>Status</b></td><td><b>Input file</b></td>
                   <td><b>Output file</b></td></tr></thead>
          <tbody>
            <xsl:for-each select="TASK">
              <tr>
                <td><xsl:value-of select="@status"/></td>
                <td>
                  <xsl:for-each select="INPUT">
                    <A href="{@file}"><xsl:value-of select="@file"/></A>
                  </xsl:for-each>
                </td>
                <td>
                  <xsl:for-each select="OUTPUT">
                    <A href="{@file}"><xsl:value-of select="@file"/> </A>
                  </xsl:for-each>
                </td>
              </tr>
            </xsl:for-each>
          </tbody>
        </table>
      </body>
    </html>
  </xsl:template>

  <xsl:template match="/SIMULATION">
    <html>
      <head/>
      <body>
        <H1>ALPS Simulation</H1>
          <H2>Parameters</H2>
          <table border="1">
            <thead><tr><td>
                    <B>Parameter</B>
                  </td>
                  <td>
                    <B>Value</B>
                  </td>
                </tr>
              </thead>
          <tbody>
          <xsl:for-each select="PARAMETERS">
            <xsl:for-each select="PARAMETER">
              <tr>
                <td>
                  <xsl:value-of select="@name"/>
                </td>
                <td>
                  <xsl:apply-templates/>
                </td>
              </tr>
            </xsl:for-each>
          </xsl:for-each>
          </tbody></table>
          <H2>Execution information</H2>
          <xsl:for-each select="MCRUN">
            <H3>Run <xsl:value-of select="position()"/></H3>
            <table border="1"><thead><tr><td><b>Executed from</b></td><td><b>to</b></td>
              <td><b>on</b></td></tr></thead><tbody>
              <xsl:for-each select="EXECUTED">
              <tr>
                <td><xsl:for-each select="FROM"><xsl:apply-templates/></xsl:for-each></td>
                <td><xsl:for-each select="TO"><xsl:apply-templates/></xsl:for-each></td>
                <td><xsl:for-each select="MACHINE">
                      <xsl:for-each select="NAME"><xsl:apply-templates/></xsl:for-each>
                    </xsl:for-each></td>
              </tr>
              </xsl:for-each>
            </tbody></table>
          </xsl:for-each>
        <H2>Averages</H2>
          <H3>Total</H3>
            <UL>
            <xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
              <LI><A HREF="#scalar{position()}"><xsl:value-of select="@name"/></A></LI>
            </xsl:for-each>
            <xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
              <LI><A HREF="#vector{position()}"><xsl:value-of select="@name"/></A></LI>
            </xsl:for-each>
            </UL>
            <TABLE BORDER="1" COLS="5" WIDTH="800">
            
            <THEAD><TR><TD><B>Name</B></TD><TD><B>Count</B></TD><TD><B>Mean</B></TD>
                       <TD><B>Error</B></TD><TD><B>Tau</B></TD>
                       <TD><B>Method</B></TD></TR></THEAD>
            <TBODY>
            <xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
            <A NAME="scalar{position()}"></A>
            <TR>
              <TD><B><xsl:value-of select="@name"/></B></TD>
              <TD><xsl:value-of select="COUNT"/></TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <xsl:choose>
                            <xsl:when test="ERROR/@converged = 'maybe'">
              <TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@converged = 'no'">
              <TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@underflow = 'true'">
              <TD bgcolor="#00ff00"><xsl:value-of select="ERROR"/><BR/>possible underflow</TD>
                            </xsl:when>
                            <xsl:otherwise>
              <TD ><xsl:value-of select="ERROR"/></TD>
                            </xsl:otherwise>
                            </xsl:choose>
                            <TD><xsl:value-of select="AUTOCORR"/></TD>
                            <TD><xsl:value-of select="ERROR/@method"/></TD>
            </TR>
            </xsl:for-each>
            <xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
            <A NAME="vector{position()}"></A>
            <xsl:for-each select="SCALAR_AVERAGE">
            <TR>
              <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of 
              select="@indexvalue"/>]</TD>
              <TD><xsl:value-of select="COUNT"/></TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <xsl:choose>
                            <xsl:when test="ERROR/@converged = 'maybe'">
              <TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@converged = 'no'">
              <TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@underflow = 'true'">
              <TD bgcolor="#00ff00"><xsl:value-of select="ERROR"/><BR/>possible underflow</TD>
                            </xsl:when>
                            <xsl:otherwise>
              <TD ><xsl:value-of select="ERROR"/></TD>
                            </xsl:otherwise>
                            </xsl:choose>
              <TD><xsl:value-of select="AUTOCORR"/></TD>
              <TD><xsl:value-of select="ERROR/@method"/></TD>
            </TR>
            </xsl:for-each>
            </xsl:for-each>
            </TBODY>
            </TABLE>
          <xsl:for-each select="MCRUN">
            <H3>Run <xsl:value-of select="position()"/></H3>
            <TABLE BORDER="1" COLS="5" WIDTH="800">
            
            <THEAD><TR><TD><B>Name</B></TD><TD><B>Count</B></TD><TD><B>Mean</B></TD>
                       <TD><B>Error</B></TD><TD><B>Tau</B></TD>
                       <TD><B>Method</B></TD></TR></THEAD>
            <TBODY>
            <xsl:for-each select="AVERAGES/SCALAR_AVERAGE">
            <TR>
              <TD><B><xsl:value-of select="@name"/></B></TD>
              <TD><xsl:value-of select="COUNT"/></TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <xsl:choose>
                            <xsl:when test="ERROR/@converged = 'maybe'">
              <TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@converged = 'no'">
              <TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                            </xsl:when>
                            <xsl:otherwise>
              <TD ><xsl:value-of select="ERROR"/></TD>
                            </xsl:otherwise>
                            </xsl:choose>
              <TD><xsl:value-of select="AUTOCORR"/></TD>
              <TD><xsl:value-of select="ERROR/@method"/></TD>
            </TR>
            <xsl:for-each select="BINNED">
            <TR>
              <TD><B><xsl:value-of select="@name"/></B></TD>
              <TD><xsl:value-of select="COUNT"/> bins</TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <TD ><xsl:value-of select="ERROR"/></TD>
              <TD/>
              <TD/>
            </TR>
            </xsl:for-each>
            </xsl:for-each>
            <xsl:for-each select="AVERAGES/VECTOR_AVERAGE">
            <xsl:for-each select="SCALAR_AVERAGE">
            <TR>
              <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of 
              select="@indexvalue"/>]</TD>
              <TD><xsl:value-of select="COUNT"/></TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <xsl:choose>
                            <xsl:when test="ERROR/@converged = 'maybe'">
              <TD bgcolor="#ffff00"><xsl:value-of select="ERROR"/><BR/>check convergence</TD>
                            </xsl:when>
                            <xsl:when test="ERROR/@converged = 'no'">
              <TD bgcolor="#ff0000"><blink><xsl:value-of select="ERROR"/><BR/>not converged</blink></TD>
                            </xsl:when>
                            <xsl:otherwise>
              <TD ><xsl:value-of select="ERROR"/></TD>
                            </xsl:otherwise>
                            </xsl:choose>
              <TD><xsl:value-of select="AUTOCORR"/></TD>
              <TD><xsl:value-of select="ERROR/@method"/></TD>
            </TR>
            <xsl:for-each select="BINNED">
            <TR>
              <TD><B><xsl:value-of select="@name"/></B></TD>
              <TD><xsl:value-of select="COUNT"/> bins</TD>
              <TD><xsl:value-of select="MEAN"/></TD>
                            <TD ><xsl:value-of select="ERROR"/></TD>
              <TD/>
              <TD/>
            </TR>
            </xsl:for-each>
            </xsl:for-each>
            </xsl:for-each>
            </TBODY>
            </TABLE>
</xsl:for-each>
        <H2>Eigenvalues</H2>
        <xsl:for-each select="EIGENVALUES">
          <xsl:for-each select="QUANTUMNUMBER">
            <B> Quantum number <xsl:value-of select="@name"/>=<xsl:value-of select="@value"/></B><BR/>
          </xsl:for-each>
          <PRE><xsl:value-of select="."/></PRE>
        </xsl:for-each>
        <H2>Eigen states</H2>
        <xsl:for-each select="EIGENSTATES">
          <H3>Sector</H3>
          <xsl:for-each select="QUANTUMNUMBER">
            <B> Quantum number <xsl:value-of select="@name"/>=<xsl:value-of select="@value"/></B><BR/>
          </xsl:for-each>
        <xsl:for-each select="EIGENSTATE">
          <H4>Eigenstate #<xsl:value-of select="position()"/></H4>
          <TABLE BORDER="1" COLS="5" WIDTH="800">
            <THEAD><TR><TD><B>Name</B></TD><TD><B>Expectation value</B></TD></TR></THEAD>
            <TBODY>
              <xsl:for-each select="SCALAR_AVERAGE">
              <TR>
                <TD><B><xsl:value-of select="@name"/></B></TD>
                <TD><xsl:value-of select="MEAN"/></TD>
              </TR>
              </xsl:for-each>
              <xsl:for-each select="VECTOR_AVERAGE">
                <xsl:for-each select="SCALAR_AVERAGE">
                <TR>
                  <TD><B><xsl:value-of select="../@name"/></B>[<xsl:value-of select="@indexvalue"/>]</TD>
                  <TD><xsl:value-of select="MEAN"/></TD>
                </TR>
              </xsl:for-each>
            </xsl:for-each>
          </TBODY>
        </TABLE>
      </xsl:for-each>
      </xsl:for-each>
      </body>
    </html>
  </xsl:template>

  <xsl:variable name="newline">
  <xsl:text>
  </xsl:text>
  </xsl:variable>
  
  <xsl:template match="plot">
  <html>
   <head/>
   <body>
    <h1><xsl:value-of select="@name"/></h1>
    <xsl:value-of select="yaxis/@label"/><xsl:text> versus </xsl:text>
    <xsl:value-of select="xaxis/@label"/>

    <xsl:for-each select="set">
      <xsl:if test= "../legend/@show = 'true' and @label">
        <h2><xsl:value-of select="@label"/></h2>
        <xsl:value-of select="$newline"/>
      </xsl:if>

      <table border="1">
      <xsl:apply-templates select="point">
        <xsl:sort select="x" data-type="number"/>
      </xsl:apply-templates>
      </table>
      
    </xsl:for-each>

   </body>
  </html>
</xsl:template>

<xsl:template match="plot/set/point">
 <xsl:if test="position()=1">
 <tr>
   <td><b><xsl:value-of select="../../xaxis/@label"/></b></td>
   <xsl:if test= "dx">
      <td><b><xsl:value-of select="../../xaxis/@label"/>
          <xsl:text> Error</xsl:text></b></td>
   </xsl:if>
   <td><b><xsl:value-of select="../../yaxis/@label"/></b></td>
   <xsl:if test= "dy">
      <td><b><xsl:value-of select="../../yaxis/@label"/>
          <xsl:text> Error</xsl:text></b></td>
   </xsl:if>
 </tr>
 </xsl:if>
 <xsl:if test="x and y">
   <tr>
     <td><xsl:value-of select="x"/></td>
     <xsl:if test= "dx">
       <td><xsl:value-of select="dx"/></td>
     </xsl:if>
     <td><xsl:value-of select="y"/></td>
     <xsl:if test= "dy">
       <td><xsl:value-of select="dy"/></td>
     </xsl:if>
   </tr> 
 </xsl:if>
</xsl:template>



</xsl:stylesheet>
