subroutine corr_coef_point(fi, prdct, polyVar, lon_size, lat_size, z_size, t_size, fi_size, corrx)
implicit none

! size of the arrays
integer, intent(in) :: lon_size,lat_size,z_size,t_size,fi_size,fi

! data to be used for the correlation calculation
real, intent(in), dimension(z_size, t_size, lat_size, lon_size) :: prdct
real, intent(in), dimension(t_size, fi_size) :: polyVar

! array to hold the outcome of the correlation
real, intent(out), dimension(z_size, lat_size, lon_size) :: corrx
!f2py intent(out) :: corrx

! variable needed for the loop
integer :: loni,lati,zi
real :: mu_polyVar,mu_prdct,covariance,std_polyVar,std_prdct,front_var,one_var

mu_polyVar = sum(polyVar(:,fi)) / t_size
one_var = 1
front_var = (one_var /(t_size - one_var))

! loooooooooop
do loni = 1, lon_size
    do lati = 1, lat_size
        do zi = 1, z_size
            mu_prdct   = sum(prdct(zi, :, lati, loni)) / t_size

            covariance = front_var * (sum((polyVar(:,fi) - mu_polyVar) * &
            (prdct(zi, :, lati, loni) - mu_prdct)))

            std_polyVar = front_var * (sum((polyVar(:,fi) - mu_polyVar) * &
            (polyVar(:,fi) - mu_polyVar)))

            std_prdct   = front_var * (sum((prdct(zi, :, lati, loni) - mu_prdct) * &
            (prdct(zi, :, lati, loni) - mu_prdct)))

            corrx(zi, lati, loni) = covariance / (SqRt(std_polyVar) * SqRt(std_prdct))

        end do
    end do
end do


end subroutine corr_coef_point