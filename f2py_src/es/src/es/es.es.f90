SUBROUTINE calculate_es()
!###########################################################
INTEGER(KIND=StandardInteger) :: calc_id
!###########################################################

! Reset total rss over all 
rss_total_rss = 0.0D0
rss_total_rss_w = 0.0D0

! Loop through
calc_id = 1
DO WHILE(calc_keys_i(calc_id, 1) .NE. -1)
  !CALL calculate_bp_inner(bp_id)
  print *, calc_id
  calc_id = calc_id + 1
END DO

END SUBROUTINE calculate_es
