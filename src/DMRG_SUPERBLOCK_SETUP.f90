!This module loads all the public procedures to setup the SuperBlock Hamiltonian
MODULE DMRG_SUPERBLOCK_SETUP
  USE DMRG_SUPERBLOCK_SETUP_GLOBAL, only: Free_SuperBlock_setup_state
  USE DMRG_SUPERBLOCK_SETUP_SPARSE, only: Setup_SuperBlock_Sparse
  USE DMRG_SUPERBLOCK_SETUP_DIRECT, only: Setup_SuperBlock_Direct
  implicit none
  private

  public :: Setup_SuperBlock_Sparse
  public :: Setup_SuperBlock_Direct
  public :: free_superblock_setup_state

END MODULE DMRG_SUPERBLOCK_SETUP
