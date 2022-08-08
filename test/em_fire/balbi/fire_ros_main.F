program fire_ros_main

! call fire_ros similarly as in write_fuels_m

! include fuel categories etc. fire_ros meeds, as module variables 
use module_fr_sfire_phys

type(fire_params)::fp  ! this is not a nodule variable we could include

! fake bounds, all is one fire mesh cell only
integer::ifds=1,ifde=1,jfds=1,jfde=1,ifms=1,ifme=1,jfms=1,jfme=1,ifts=1,ifte=1,jfts=1,jfte=1

! read namelist.fire, do not initialize the fuel moisture model
call init_fuel_cats(.false.) 

! read the input variables as namelist

call set_fire_params( &
                           ifds,ifde,jfds,jfde, &
                           ifms,ifme,jfms,jfme, &
                           ifts,ifte,jfts,jfte, &
                           fdx,fdy,nfuel_cat0,  &
                           nfuel_cat,fuel_time, &
                           fp )


end module module_call_fire_ros




end program fire_ros_main
