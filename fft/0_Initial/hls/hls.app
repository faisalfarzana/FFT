<AutoPilot:project xmlns:AutoPilot="com.autoesl.autopilot.project" projectType="C/C++" name="hls" ideType="classic" top="fft">
    <files>
        <file name="fft.h" sc="0" tb="false" cflags="" csimflags="" blackbox="false"/>
        <file name="fft.cpp" sc="0" tb="false" cflags="" csimflags="" blackbox="false"/>
        <file name="../../out.gold.dat" sc="0" tb="1" cflags="-Wno-unknown-pragmas" csimflags="" blackbox="false"/>
        <file name="../../fft_test.cpp" sc="0" tb="1" cflags="-Wno-unknown-pragmas" csimflags="" blackbox="false"/>
    </files>
    <solutions>
        <solution name="solution1" status=""/>
    </solutions>
    <Simulation argv="">
        <SimFlow name="csim" setup="false" optimizeCompile="false" clean="true" ldflags="" mflags=""/>
    </Simulation>
</AutoPilot:project>

