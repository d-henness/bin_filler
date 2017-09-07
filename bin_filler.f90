module find_bin_configs
  implicit none
  private
  public  ::  find_bin_configurations
  contains
    subroutine find_bin_configurations(nbins, nactive_bins, active_bin_list, bin_sizes, nobjects, nconfigs, all_configs)
    ! This takes any set of open and closed bins, ordered from smallest to largest, as well as a number of
    ! equally sized objects and finds all possible ways that the objects could fit
    ! into the open bins.
      integer ::  nbins, nactive_bins, nobjects, i, j, nconfigs, istat
      integer ::  highest_occ, next_not_full, pivot, left_over
      integer, dimension(:) ::  active_bin_list, bin_sizes
      integer, dimension(:), allocatable  ::  bin_status,  bin_mapping
      real(kind=8), dimension(:, :), allocatable ::  all_configs
      
      allocate(bin_status(nbins), stat = istat)
      if (istat /= 0) stop 'Error during allocation of bin_status'
      allocate(bin_mapping(nactive_bins), stat = istat)
      if (istat /= 0) stop 'Error during allocation of bin_mapping'
    
      j = 1
      do i = 1, nbins
        if (active_bin_list(i) == 1) then
          bin_mapping(j) = i
          j = j + 1
        endif
      enddo
    
      if (j - 1 /= nactive_bins) stop 'Number of active bins does not match the active_bin_list'
    
      bin_status = 0
      nconfigs = 0

      ! I can't think of a way to figure out how many configurations there will
      ! be without just finding them all first. If you can, just replace from
      ! here until the stop comment

      call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, nobjects, nconfigs)
    
      do while(.true.)
        highest_occ = 1
        next_not_full = 0
        ! Find the largest occcupied bin
        do i = nactive_bins, 1, -1
          if (bin_status(bin_mapping(i)) > 0) then 
            highest_occ = i
            exit
          endif
        enddo
        if (highest_occ == 1) exit ! If the highest occupied bin is the smallest bin, then all configurations have been found
        ! Find the next bin that isn't full
        do i = highest_occ - 1, 1, -1
          if (bin_status(bin_mapping(i)) < bin_sizes(bin_mapping(i))) then
            next_not_full = i
            exit
          endif
        enddo
        if (next_not_full == 0) then 
          exit ! If there are not any bins smaller than the highest occupied bin that are also not full, then all configs have been found
        else if (next_not_full == 1) then ! If this is the case, we add one to the next_not_full bin and then reset the config minus this object
          bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) + 1
          bin_status(bin_mapping(highest_occ)) = bin_status(bin_mapping(highest_occ)) - 1
          left_over = 0
          do i = next_not_full + 1, nactive_bins
            left_over = left_over + bin_status(bin_mapping(i))
            bin_status(bin_mapping(i)) = 0
          enddo
          call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, left_over, nconfigs)
          cycle
        endif
        ! Move as many objects from highest_occ bin to the next_not_full bin as are allowed
        if (next_not_full + 1 == highest_occ) then ! If the next_not_full and highest_occ are next to each other, move objects
          do while ((bin_status(bin_mapping(highest_occ)) /= 0) .and. (bin_status(bin_mapping(next_not_full)) < bin_sizes(bin_mapping(next_not_full))))
            bin_status(bin_mapping(highest_occ)) = bin_status(bin_mapping(highest_occ)) - 1
            bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) + 1
            nconfigs = nconfigs + 1
          enddo
        else ! Otherwise, we need to change the next_not_full to the highest_occ for the pivot finder to work
          next_not_full = highest_occ
        endif
        ! Now we have to find the next_not_full bin (pivot) after the current next_not_full bin.
        pivot = 0
        do i = next_not_full - 1, 1, -1
          if (bin_status(bin_mapping(i)) < bin_sizes(bin_mapping(i))) then
            pivot = i
            exit
          endif
        enddo
        if (pivot == 0) exit ! If there are not any bins smaller than the next_not_full bin that are also not full, then all configs have been found
        ! Add one to the pivot, find out how many objects are left, and then resets the bin occupation minus the objects in the pivot bin
        bin_status(bin_mapping(pivot)) = bin_status(bin_mapping(pivot)) + 1
        bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) - 1
        left_over = 0
        do i = pivot + 1, nactive_bins
          left_over = left_over + bin_status(bin_mapping(i))
          bin_status(bin_mapping(i)) = 0
        enddo
        call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, left_over, nconfigs)
      enddo

      !!!!!!!!!!!!!!!!!!!!!
      ! Replace till here !
      !!!!!!!!!!!!!!!!!!!!!

      write(6, *)'nconfigs is ', nconfigs

      allocate(all_configs(nbins, nconfigs), stat = istat)
      if (istat /= 0) stop 'Error during allocation of all_configs'

      nconfigs = 0
      bin_status = 0
      call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, nobjects, nconfigs)
      all_configs(:, nconfigs) = bin_status
    
      do while(.true.)
        highest_occ = 1
        next_not_full = 0
        ! Find the largest occcupied bin
        do i = nactive_bins, 1, -1
          if (bin_status(bin_mapping(i)) > 0) then 
            highest_occ = i
            exit
          endif
        enddo
        if (highest_occ == 1) exit ! If the highest occupied bin is the smallest bin, then all configurations have been found
        ! Find the next bin that isn't full
        do i = highest_occ - 1, 1, -1
          if (bin_status(bin_mapping(i)) < bin_sizes(bin_mapping(i))) then
            next_not_full = i
            exit
          endif
        enddo
        if (next_not_full == 0) then 
          exit ! If there are not any bins smaller than the highest occupied bin that are also not full, then all configs have been found
        else if (next_not_full == 1) then ! If this is the case, we add one to the next_not_full bin and then reset the config minus this object
          bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) + 1
          bin_status(bin_mapping(highest_occ)) = bin_status(bin_mapping(highest_occ)) - 1
          left_over = 0
          do i = next_not_full + 1, nactive_bins
            left_over = left_over + bin_status(bin_mapping(i))
            bin_status(bin_mapping(i)) = 0
          enddo
          call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, left_over, nconfigs)
          all_configs(:, nconfigs) = bin_status
          cycle
        endif
        ! Move as many objects from highest_occ bin to the next_not_full bin as are allowed
        if (next_not_full + 1 == highest_occ) then ! If the next_not_full and highest_occ are next to each other, move objects
          do while ((bin_status(bin_mapping(highest_occ)) /= 0) .and. (bin_status(bin_mapping(next_not_full)) < bin_sizes(bin_mapping(next_not_full))))
            bin_status(bin_mapping(highest_occ)) = bin_status(bin_mapping(highest_occ)) - 1
            bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) + 1
            nconfigs = nconfigs + 1
            all_configs(:, nconfigs) = bin_status
          enddo
        else ! Otherwise, we need to change the next_not_full to the highest_occ for the pivot finder to work
          next_not_full = highest_occ
        endif
        ! Now we have to find the next_not_full bin (pivot) after the current next_not_full bin.
        pivot = 0
        do i = next_not_full - 1, 1, -1
          if (bin_status(bin_mapping(i)) < bin_sizes(bin_mapping(i))) then
            pivot = i
            exit
          endif
        enddo
        if (pivot == 0) exit ! If there are not any bins smaller than the next_not_full bin that are also not full, then all configs have been found
        ! Add one to the pivot, find out how many objects are left, and then resets the bin occupation minus the objects in the pivot bin
        bin_status(bin_mapping(pivot)) = bin_status(bin_mapping(pivot)) + 1
        bin_status(bin_mapping(next_not_full)) = bin_status(bin_mapping(next_not_full)) - 1
        left_over = 0
        do i = pivot + 1, nactive_bins
          left_over = left_over + bin_status(bin_mapping(i))
          bin_status(bin_mapping(i)) = 0
        enddo
        call fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, left_over, nconfigs)
        all_configs(:, nconfigs) = bin_status
      enddo
    end subroutine find_bin_configurations

    subroutine fill_bins(bin_status, bin_sizes, bin_mapping, nactive_bins, nobjects, nconfigs)
      integer, dimension(:) ::  bin_status, bin_sizes, bin_mapping
      integer ::  nobjects, nactive_bins, nconfigs
      integer ::  objects_left, indx
    
      objects_left = nobjects
      indx = nactive_bins
      do while ((objects_left /= 0) .and. (indx /= 0))
        if (bin_status(bin_mapping(indx)) < bin_sizes(bin_mapping(indx))) then
          bin_status(bin_mapping(indx)) = bin_status(bin_mapping(indx)) + 1
          objects_left = objects_left - 1
        else
          indx = indx - 1
        endif
      enddo
    
      if (indx == 0) stop 'Too many objects for number and size of bins'
      nconfigs = nconfigs + 1
    end subroutine fill_bins
end module find_bin_configs

program main
  use find_bin_configs
  implicit none

  integer ::  nbins, nactive_bins, nobjects, i, j, nconfigs, istat
  integer, allocatable, dimension(:) ::  active_bin_list, bin_sizes
  real(kind=8), dimension(:, :), allocatable ::  all_configs

  read(5, *)nbins, nactive_bins

  allocate(bin_sizes(nbins), active_bin_list(nbins))

  read(5, *)active_bin_list
  read(5, *)bin_sizes
  read(5, *)nobjects

  call find_bin_configurations(nbins, nactive_bins, active_bin_list, bin_sizes, nobjects, nconfigs, all_configs)

  write(6, *)'# of configis', nconfigs
  do i = 1, nconfigs
    write(6, fmt='(7(f3.0))')all_configs(:, i)
  enddo

end program main
