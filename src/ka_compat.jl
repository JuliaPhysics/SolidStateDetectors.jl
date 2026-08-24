# This file is a part of RadiationDetectorDSP.jl, licensed under the MIT License (MIT).


@static if isdefined(KernelAbstractions, :Backend)
    import KernelAbstractions.Backend as _KA_Backend
    import KernelAbstractions.get_backend as _ka_get_backend
    # GPU backends are stream-ordered: subsequent operations on the same queue see
    # the kernel's results without a host-side synchronization. Only the CPU
    # backend (task-based) needs an explicit wait.
    _ka_synchronize(kernel, kernel_ret) = kernel.backend isa KernelAbstractions.CPU ?
        KernelAbstractions.synchronize(kernel.backend) : nothing
else
    # KernelAbstractions < v0.9:
    import KernelAbstractions.Device as _KA_Backend
    import KernelAbstractions.get_device as _ka_get_backend
    _ka_synchronize(kernel, kernel_ret) = wait(kernel_ret)
end
