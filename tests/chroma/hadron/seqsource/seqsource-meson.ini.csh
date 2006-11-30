#!/bin/csh

cat <<EOF
<?xml version="1.0"?>
<chroma>
<annotation>
Sequential source
</annotation>
<Param> 
  <InlineMeasurements>

    <elem>
      <Name>MAKE_SOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>SHELL_SOURCE</SourceType>
          <j_decay>3</j_decay>
          <t_srce>0 0 0 6</t_srce>

          <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>

          <LinkSmearing>
            <LinkSmearingType>APE_SMEAR</LinkSmearingType>
            <link_smear_fact>2.5</link_smear_fact>
            <link_smear_num>0</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Source>

      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_0</source_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
         <FermAct>WILSON</FermAct>
         <Kappa>0.11</Kappa>
         <AnisoParam>
           <anisoP>false</anisoP>
           <t_dir>3</t_dir>
           <xi_0>1.0</xi_0>
           <nu>1.0</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>1 1 1 0</boundary>
         </FermionBC>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>1.0e-12</RsdCG>
          <MaxCG>1000</MaxCG>
        </InvertParam>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>sh_source_0</source_id>
        <prop_id>sh_prop_0</prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SINK_SMEAR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>5</version>
        <Sink>
          <version>1</version>
          <SinkType>SHELL_SINK</SinkType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>

          <LinkSmearing>
            <LinkSmearingType>APE_SMEAR</LinkSmearingType>
            <link_smear_fact>2.5</link_smear_fact>
            <link_smear_num>0</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Sink>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_id>sh_prop_0</prop_id>
        <smeared_prop_id>sh_pt_prop_0</smeared_prop_id>
      </NamedObject>
    </elem>

EOF

foreach seq_src (PION_1-PION_1 A0-A0_1 A0-RHO_X_1 A0-RHO_Y_1 A0-B1_Z_1 A0-RHO_Z_1 A0-B1_Y_1 A0-B1_X_1 A0-PION_2 A0-A0_2 A0-RHO_X_2 A0-RHO_Y_2 A0-A1_Z_1 A0-RHO_Z_2 A0-A1_Y_1 A0-A1_X_1 A0-PION_1)

cat <<EOF
    <elem>
      <annotation>
       ${seq_src} seqsource
      </annotation>

      <Name>SEQSOURCE</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>1</version>
        <seq_src>${seq_src}</seq_src>
        <t_sink>6</t_sink>
        <sink_mom>1 0 0</sink_mom>
      </Param>
      <PropSink>
        <version>5</version>
        <Sink>
          <version>2</version>
          <SinkType>SHELL_SINK</SinkType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>

          <LinkSmearing>
            <LinkSmearingType>APE_SMEAR</LinkSmearingType>
            <link_smear_fact>2.5</link_smear_fact>
            <link_smear_num>0</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Sink>
      </PropSink>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <prop_ids>
          <elem>sh_prop_0</elem>
        </prop_ids>
        <seqsource_id>seqsource_tmp</seqsource_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>PROPAGATOR</Name>
      <Frequency>1</Frequency>
      <Param>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
         <FermAct>WILSON</FermAct>
         <Kappa>0.11</Kappa>
         <AnisoParam>
           <anisoP>false</anisoP>
           <t_dir>3</t_dir>
           <xi_0>1.0</xi_0>
           <nu>1.0</nu>
         </AnisoParam>
         <FermionBC>
           <FermBC>SIMPLE_FERMBC</FermBC>
           <boundary>1 1 1 0</boundary>
         </FermionBC>
        </FermionAction>
        <InvertParam>
          <invType>CG_INVERTER</invType>
          <RsdCG>1.0e-12</RsdCG>
          <MaxCG>1000</MaxCG>
        </InvertParam>
      </Param>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <source_id>seqsource_tmp</source_id>
        <prop_id>seqprop_0</prop_id>
      </NamedObject>
    </elem>

    <elem>
      <Name>SEQPROP_TEST</Name>
      <Frequency>1</Frequency>
      <PropSourceSmear>
        <version>6</version>
        <Source>
          <version>2</version>
          <SourceType>SHELL_SOURCE</SourceType>
          <j_decay>3</j_decay>

          <SmearingParam>
            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>
            <wvf_param>2.0</wvf_param>
            <wvfIntPar>5</wvfIntPar>
            <no_smear_dir>3</no_smear_dir>
          </SmearingParam>

          <Displacement>
            <version>1</version>
            <DisplacementType>NONE</DisplacementType>
          </Displacement>

          <LinkSmearing>
            <LinkSmearingType>APE_SMEAR</LinkSmearingType>
            <link_smear_fact>2.5</link_smear_fact>
            <link_smear_num>0</link_smear_num>
            <no_smear_dir>3</no_smear_dir>
          </LinkSmearing>
        </Source>
      </PropSourceSmear>
      <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <sink_ids>
          <elem>sh_pt_prop_0</elem>
          <elem>sh_pt_prop_0</elem>
        </sink_ids>
        <seqprop_id>seqprop_0</seqprop_id>
        <gamma_insertion>0</gamma_insertion>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
        Erase the object
      </annotation>
      <Name>ERASE_NAMED_OBJECT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <object_id>seqsource_tmp</object_id>
      </NamedObject>
    </elem>

    <elem>
      <annotation>
        Erase the object
      </annotation>
      <Name>ERASE_NAMED_OBJECT</Name>
      <Frequency>1</Frequency>
      <NamedObject>
        <object_id>seqprop_0</object_id>
      </NamedObject>
    </elem>

EOF
end  # foreach

cat <<EOF
  </InlineMeasurements>
   <nrow>4 4 4 8</nrow>
</Param>

<RNG>
  <Seed>	
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

<Cfg>
 <cfg_type>WEAK_FIELD</cfg_type>
 <cfg_file>dummy</cfg_file>
</Cfg>
</chroma>

EOF

