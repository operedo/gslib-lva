      module kind
        !INTEGER, PARAMETER :: DP=KIND(1.0D0)
        integer, parameter :: sp = selected_real_kind(6, 37)
        !integer, parameter :: dp = selected_real_kind(6, 37)
        integer, parameter :: dp = selected_real_kind(15, 307)
        !integer, parameter :: qp = selected_real_kind(33, 4931)
      end module kind
