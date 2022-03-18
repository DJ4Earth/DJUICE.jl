function cuffey(temperature)
#CUFFEY - calculates ice rigidity as a function of temperature
#
#   rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3)
#   (Cuffey and Paterson, p75).
#   temperature is in Kelvin degrees

#   Usage:
#      rigidity=cuffey(temperature)

	if temperature<0
	   error("input temperature should be in Kelvin (positive)")
	end
	T=temperature-273.15

#The routine below is equivalent to:
# n=3; T=temperature-273;
# %From cuffey
# Temp=[0;-2;-5;-10;-15;-20;-25;-30;-35;-40;-45;-50];
# A=[2.4*10^-24;1.7*10^-24;9.3*10^-25;3.5*10^-25;2.1*10^-25;1.2*10^-25;6.8*10^-26;3.7*10^-26;2.0*10^-26;1.0*10^-26;5.2*10^-27;2.6*10^-27];%s-1(Pa-3)
# %Convert into rigidity B
# B=A.^(-1/n); %s^(1/3)Pa
# %Now, do a cubic fit between Temp and B:
# fittedmodel=fit(Temp,B,'cubicspline');
# rigidity=fittedmodel(temperature);

	rigidity=zeros(length(T));

	if T<=-45
		rigidity = 10^8*(-0.000396645116301*(T+50).^3+ 0.013345579471334*(T+50).^2  -0.356868703259105*(T+50)+7.272363035371383)
	elseif -45<=T && T<-40
		rigidity=10^8*(-0.000396645116301*(T+45).^3+ 0.007395902726819*(T+45).^2  -0.253161292268336*(T+45)+5.772078366321591)
	elseif -45<=T && T<-40
		rigidity=10^8*(-0.000396645116301*(T+45).^3+ 0.007395902726819*(T+45).^2  -0.253161292268336*(T+45)+5.772078366321591)
	elseif -40<=T && T<-35
		rigidity=10^8*(0.000408322072669*(T+40).^3+  0.001446225982305*(T+40).^2  -0.208950648722716*(T+40)+4.641588833612773)
	elseif -35<=T && T<-30
		rigidity=10^8*(-0.000423888728124*(T+35).^3+ 0.007571057072334*(T+35).^2  -0.163864233449525*(T+35)+3.684031498640382)
	elseif -30<=T && T<-25
		rigidity=10^8*(0.000147154327025*(T+30).^3+ 0.001212726150476*(T+30).^2  -0.119945317335478*(T+30)+3.001000667185614)
	elseif -25<=T && T<-20
		rigidity=10^8*(-0.000193435838672*(T+25).^3+ 0.003420041055847*(T+25).^2  -0.096781481303861*(T+25)+2.449986525148220)
	elseif -20<=T && T<-15
		rigidity=10^8*(0.000219771255067*(T+20).^3+  0.000518503475772*(T+20).^2  -0.077088758645767*(T+20)+2.027400665191131)
	elseif -15<=T && T<-10
		rigidity=10^8*(-0.000653438900191*(T+15).^3+ 0.003815072301777*(T+15).^2  -0.055420879758021*(T+15)+1.682390865739973)
	elseif -10<=T && T<-5
		rigidity=10^8*(0.000692439419762*(T+10).^3 -0.005986511201093 *(T+10).^2 -0.066278074254598*(T+10)+1.418983411970382)
	elseif -5<=T && T<-2
		rigidity=10^8*(-0.000132282004110*(T+5).^3 +0.004400080095332*(T+5).^2    -0.074210229783403*(T+5)+ 1.024485188140279)
	elseif -2<=T
		rigidity=10^8*(-0.000132282004110*(T+2).^3 +0.003209542058346*(T+2).^2    -0.051381363322371*(T+2)+ 0.837883605537096)
	end

	#Make sure that rigidity is positive
	if rigidity<0
		rigidity=10^6
	end

	return rigidity
end
