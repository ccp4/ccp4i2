/*
<CCP4i2ReportFlotGraph key="SummaryGraph" class="" style="height:250px; width:400px;float:left;border:0px;">
        <ns0:ccp4_data title="Running refmac" id="data_SummaryGraph" style="display:none;">
          <headers>Cycle R_Factor R_Free rmsBonds </headers>
          <data>0 0.2631 - 0.010 
1 0.2500 - 0.009 
2 0.2440 - 0.009 
3 0.2407 - 0.010 
4 0.2386 - 0.010 
5 0.2373 - 0.010 
6 0.2363 - 0.010 
7 0.2356 - 0.010 
8 0.2350 - 0.010 
9 0.2345 - 0.011 
10 0.2341 - 0.011 
</data>
          <plot>
            <title>Running refmac R-factors</title>
            <plottype>xy</plottype>
            <yrange rightaxis="false"/>
            <xlabel>Cycle</xlabel>
            <xintegral>true</xintegral>
            <ylabel>R-factor</ylabel>
            <rylabel>Geometry</rylabel>
            <plotline xcol="1" ycol="2" rightaxis="false" colour="blue"/>
            <plotline xcol="1" ycol="3" rightaxis="false" colour="green"/>
            <yrange rightaxis="true"/>
            <plotline xcol="1" ycol="4" rightaxis="true" colour="red"/>
          </plot>
        </ns0:ccp4_data>
      </CCP4i2ReportFlotGraph>
      */
interface ReportGraphWithDataProps {
  // Define the props for the component if needed
}
export const ReportGraphWithData: React.FC<ReportGraphWithDataProps> = ({}) => {
  return <span>Hello</span>;
};
