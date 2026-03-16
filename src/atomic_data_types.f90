! Shared atomic data types for Kurucz and NIST DB (avoids circular dependency).
      module atomic_data_types
      implicit none
      integer , parameter :: max_ion_levels=8
      integer , parameter :: max_atoms=99

      type level_data
        real(16) :: g=0.0d0
        real(16) :: energy=1.d90
      end type level_data

      type line_data
        integer :: z=0
        integer :: ion=0
        real(16) :: lambda=0.0d0
        real(16) :: elo=0.0d0
        integer :: ibin=0
        integer :: l1=0
        integer :: l2=0
        real(16) :: fij=0.0d0
      end type line_data

      type ion_data
        integer :: nlevels=0
        integer :: nlines=0
        real(8) , allocatable :: partition(:)
        real(8) , allocatable :: ptemp(:)
        type (level_data) , allocatable :: level(:)
      end type ion_data

      type atom_data
        logical :: active=.false.
        integer :: z=0
        real(8) :: ionization_energy(0:max_ion_levels)=0.0d0
        type (ion_data) :: ion(0:max_ion_levels)
      end type atom_data

      end module atomic_data_types
