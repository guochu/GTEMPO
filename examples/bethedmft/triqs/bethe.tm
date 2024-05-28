<TeXmacs|2.1>

<style|generic>

<\body>
  In the Bethe lattice case, the self consistency condition is simply (see
  <slink|https://triqs.github.io/cthyb/latest/guide/dmft.html> and
  <slink|https://triqs.github.io/triqs/latest/userguide/python/dmft_one_page.html>)

  <\equation>
    G<rsub|0><rsup|-1><around*|(|i\<omega\><rsub|n>|)>=i\<omega\><rsub|n>+\<mu\>-t<rsup|2>G<around*|(|i\<omega\><rsub|n>|)>,
  </equation>

  where <math|2t> is the half bandwidth and <math|t=0.5>. This means the
  hybridization function is simply

  <\equation>
    \<Delta\><around*|(|i\<omega\><rsub|n>|)>=t<rsup|2>G<around*|(|i\<omega\><rsub|n>|)>.
  </equation>

  The DMFT procedures are following:

  <\enumerate>
    <item>At the initial time, set <math|J<around*|(|\<varepsilon\>|)>=<frac|1|2\<pi\>t><sqrt|1-<around*|(|\<varepsilon\>/2t|)><rsup|2>>>
    and we have

    <\equation>
      \<Delta\><around*|(|i\<omega\><rsub|n>|)>=<big|int><rsub|-2t><rsup|2t><frac|J<around*|(|\<varepsilon\>|)>|i\<omega\><rsub|n>-\<varepsilon\>>\<mathd\>\<varepsilon\>.
    </equation>

    <item>Compute <math|G<around*|(|i\<omega\><rsub|n>|)>>, and set new
    <math|\<Delta\><around*|(|i\<omega\><rsub|n>|)>=t<rsup|2>G<around*|(|i\<omega\><rsub|n>|)>>
  </enumerate>

  The number of mesh points for <math|i\<omega\><rsub|s>> is totally 2050 for
  both positive and negative parts, and number of mesh points for
  <math|\<tau\>> is 10001.
</body>

<\initial>
  <\collection>
    <associate|font|math=termes,times>
    <associate|font-family|rm>
  </collection>
</initial>