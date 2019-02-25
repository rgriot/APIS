subroutine likelihoodCalculation(offGeno, sires, dams, nMrk, nVariant, nSires, nDams, F, &
																 output_sires, output_dams, output_score, output_miss)

	!DEC$ ATTRIBUTES DLLEXPORT :: likelihoodCalculation
	integer			 															:: nMrk, nSires, nDams, sID, dID
	integer            												:: i, j, k, line
	double precision, parameter    						:: e = 0.01
  double precision, dimension(nMrk)  				:: likelihood, f0, f1, logLikelihood
	double precision                      		:: l, fa, fb, fc
	double precision, dimension(16)						:: homoScore, homoMiss
	double precision, dimension(7*7)					:: heteroScore, heteroMiss
	double precision, dimension(4,4) 					:: tableHomo, table_homoMiss
	double precision, dimension(7,7) 					:: tableHetero, table_heteroMiss
	double precision, dimension(nMrk,nVariant):: F
  integer, dimension(nMrk)          				:: mrkPosition, mismatch

	integer, dimension(2*nMrk) 			  				:: offGeno

	integer, dimension(nSires,2*nMrk) 				:: sires
	integer, dimension(nDams,2*nMrk)  				:: dams

	integer, dimension(2)											:: tmpOff, tmpSire, tmpDam

	integer, dimension(nSires*nDams) 					:: output_sires, output_dams
	double precision, dimension(nSires*nDams) :: output_score
	integer, dimension(nSires*nDams)					:: output_miss

	integer 																	:: r

  mrkPosition = (/(r, r=1, 2*nMrk, 2)/)
	line = 1

	do i = 1,nSires
		do j = 1,nDams
				do k = 1,nMrk
					tmpOff 	= (/ offGeno(mrkPosition(k)), offGeno(mrkPosition(k)+1) /)
					tmpSire = (/ sires(i, mrkPosition(k)), sires(i, mrkPosition(k)+1) /)
					tmpDam  = (/ dams(j, mrkPosition(k)), dams(j, mrkPosition(k)+1) /)

					if (tmpOff(1) == 0 .and. tmpOff(2) == 0) THEN 					! If the offspring has no genotype
						likelihood(k) = 1.0d0
						mismatch(k) 	= 0
					else if (tmpOff(1) == tmpOff(2) .and. tmpOff(1) /= 0) THEN		! If the offspring is homozygous
						if (tmpSire(1) == tmpSire(2)) THEN								! If the SIRE is homozygous
							if (tmpSire(1) == 0) THEN
								sID = 4
							else if (tmpSire(1) == tmpOff(1)) THEN
								sID = 1
							else
								sID = 3
							end if
						else
							if ((tmpSire(1) /= tmpOff(1) .and. tmpSire(2) /= tmpOff(1)) &
							.and. (tmpSire(1) /= tmpOff(2) .and. tmpSire(2) /= tmpOff(2))) THEN
								sID = 3
							else
								sID = 2
							end if
						end if

						if (tmpDam(1) == tmpDam(2)) THEN								! If the Dam is homozygous
							if (tmpDam(1) == 0) THEN
								dID = 4
							else if (tmpDam(1) == tmpOff(1)) THEN
								dID = 1
							else
								dID = 3
							end if
						else
							if ((tmpDam(1) /= tmpOff(1) .and. tmpDam(2) /= tmpOff(1)) &
							.and. (tmpDam(1) /= tmpOff(2) .and. tmpDam(2) /= tmpOff(2))) THEN
								dID = 3
							else
								dID = 2
							end if
						end if
					else
						if (tmpSire(1) == tmpSire(2)) THEN
							if (tmpSire(1) == 0) THEN
								sID = 7
							else if (tmpSire(1) == tmpOff(1)) THEN
								sID = 1
							else if (tmpSire(1) == tmpOff(2)) THEN
								sID = 3
							else
								sID = 6
							end if
						else
							if ((tmpSire(1) == tmpOff(1) .or. tmpSire(2) == tmpOff(1)) &
							.and. (tmpSire(1) == tmpOff(2) .or. tmpSire(2) == tmpOff(2))) THEN
							  sID = 2
							else if ((tmpSire(1) == tmpOff(1) .or. tmpSire(2) == tmpOff(1)) &
							.and. (tmpSire(1) /= tmpOff(2) .or. tmpSire(2) /= tmpOff(2))) THEN
							  sID = 4
							else if ((tmpSire(1) /= tmpOff(1) .or. tmpSire(2) /= tmpOff(1)) &
							.and. (tmpSire(1) == tmpOff(2) .or. tmpSire(2) == tmpOff(2))) THEN
							  sID = 5
							else
							  sID = 6
							end if
						end if

						if (tmpDam(1) == tmpDam(2)) THEN
							if (tmpDam(1) == 0) THEN
								dID = 7
							else if (tmpDam(1) == tmpOff(1)) THEN
								dID = 1
							else if (tmpDam(1) == tmpOff(2)) THEN
								dID = 3
							else
								dID = 6
							end if
						else
							if ((tmpDam(1) == tmpOff(1) .or. tmpDam(2) == tmpOff(1)) &
							.and. (tmpDam(1) == tmpOff(2) .or. tmpDam(2) == tmpOff(2))) THEN
							  dID = 2
							else if ((tmpDam(1) == tmpOff(1) .or. tmpDam(2) == tmpOff(1)) &
							.and. (tmpDam(1) /= tmpOff(2) .or. tmpDam(2) /= tmpOff(2))) THEN
							  dID = 4
							else if ((tmpDam(1) /= tmpOff(1) .or. tmpDam(2) /= tmpOff(1)) &
							.and. (tmpDam(1) == tmpOff(2) .or. tmpDam(2) == tmpOff(2))) THEN
							  dID = 5
							else
							  dID = 6
							end if
						end if
					end if

					if ((tmpOff(1) == tmpOff(2)) .and. tmpOff(1) /= 0) THEN
						fa = F(k, tmpOff(1))
						homoScore = (/ 1.0d0, 0.5d0, e, fa ,0.5d0 ,0.25d0 ,e ,0.5d0*fa ,&
												   e, e, e, e, fa, 0.5d0*fa, e, fa*fa /)
						tableHomo = reshape(homoScore, (/4,4/))

						likelihood(k) = tableHomo(dID, sID)

						homoMiss = (/ 0, 0, 1, 0 ,0 ,0 ,1 ,0 ,&
												  1, 1, 2, 1, 0, 0, 1, 0 /)
						table_homoMiss = reshape(homoMiss, (/4,4/))

						mismatch(k) = table_homoMiss(dID, sID)

					else if ((tmpOff(1) == tmpOff(2)) .and. tmpOff(1) == 0) THEN
						likelihood(k) = 1.0d0
						mismatch(k) 	= 0

					else
							fa = F(k, tmpOff(1))
							fb = F(k, tmpOff(2))
						heteroScore = (/ e, 0.5d0, 1.0d0, e, 0.5d0, e, fb, 0.5d0, 0.5d0, 0.5d0, 0.25d0, 0.5d0, e, 0.5d0*(fa+fb),&
                             1.0d0, 0.5d0, e, 0.5d0, e, e, fa, e, 0.25d0, 0.5d0, e, 0.25d0, e, 0.5d0*fb,&
                             0.5d0, 0.25d0, e, 0.25d0, e, e, 0.5d0*fb, e, e, e, e, e, e, e,&
                             fb,0.5d0*(fa+fb),fa,0.5d0*fb,0.5d0*fa,e,2.0d0*fa*fb /)
						tableHetero = reshape(heteroScore, (/7,7/))

            likelihood(k) = tableHetero(dID, sID)

						heteroMiss = (/ 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,&
                            0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0,&
                            0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1,&
                            0, 0, 0, 0, 0, 1, 0 /)
						table_heteroMiss = reshape(heteroMiss, (/7,7/))

						mismatch(k) = table_heteroMiss(dID, sID)

						end if

				end do

			output_sires(line) = i
			output_dams(line)  = j

			output_score(line) = sum(log(likelihood))
			output_miss(line)	 = sum(mismatch)

			line = line + 1

		end do
	end do

end subroutine likelihoodCalculation
