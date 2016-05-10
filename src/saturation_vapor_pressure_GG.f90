
! returns the saturation vapor pressure 
! (in hPa) over water (DEPENDING ON THE TEMPERATURE) given the
! temperature(in K) %% Goff Gratch formation

subroutine saturation_vapor_pressure_GG(T,e_sat)    

real T, e_sat

integer n

   if(T>=273) then
	
	e_sat = 10 ** ( &
		 -7.90298*(373.16/T-1) &
		 + 5.02808*log10(373.16/T)  &
		 - (1.3816E-7)*(10**(11.344*(1-T/373.16)) - 1) &
		 + (8.1328E-3)*(10**(-3.49149*(373.16/T -1)) - 1) &
		 + log10(1013.246 ) )
   else
	e_sat = 10 ** ( &
		 -9.09718*(273.16/T-1) &
		 - 3.56654*log10(273.16/T) &
		 + 0.876793*(1-T/273.16) &
		 + log10(6.1071 ) );
   endif

end subroutine saturation_vapor_pressure_GG