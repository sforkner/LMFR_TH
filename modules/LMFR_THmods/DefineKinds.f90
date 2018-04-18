MODULE  DefineKinds
  !------------------------------------------------------------------
  ! Module DefineKinds Description:
  !    Define the integer,real, logical and character string types. Also
  !    define user type FileInfo to store file information.
  !    This version for GNU & Intel Fortran, PGI does not support Qkind
  !
  ! Language: Fortran 90 +
  !
  ! Author   : Samuel L. Forkner
  ! Company  : Signal Mountain Software
  ! Date     : 2018/04/16 @ 9:06 AM
  !------------------------------------------------------------------

  INTEGER (selected_int_KIND(5)), PARAMETER ::                             &
       &                              I2kind = selected_int_KIND(1)      , &
       &                              I4kind = selected_int_KIND(5)      , &
       &                              Ikind  = I4kind                    , &
       &                              L2kind = 2                         , &
       &                              L4kind = 4                         , &
       &                              Lkind  = L4Kind                    , &
       &                              Rkind = selected_real_KIND(6,37)   , &
       &                              Dkind = selected_real_KIND(15,307) , &
       &                              Qkind = selected_real_KIND(33,4931), &
       &                              Xkind = Rkind                      , &
       &                              C4kind = KIND('ABCD')              , &
       &                              C8kind = KIND('ABCDEFGH')          , &
       &                              C10kind = KIND('ABCDEFGHIJ')       , &
       &                              C12kind = KIND('ABCDEFGHIJKL')     , &
       &                              C16kind = KIND('ABCDEFGHIJKLMNOP') , &
       &                              Ckind   = C8kind
       !
       TYPE FileInfo
          CHARACTER (Len=80) name   !file name
          INTEGER (Ikind) unitno    !unit number to for file
          LOGICAL (Lkind) opened    !open status TRUE/FALSE
       END TYPE FileInfo
       !
END MODULE DefineKinds
