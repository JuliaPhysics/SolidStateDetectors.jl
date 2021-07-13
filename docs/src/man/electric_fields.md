# Electric Fields

After the electric potential is calculated, the electric field is calculated for each grid point ``(i,j,k)``:

```math
\vec{E}^{i,j,k} = \left( \mathcal{E}_r^{i,j,k}, \mathcal{E}_{\varphi}^{i,j,k}, \mathcal{E}_z^{i,j,k} \right)^{\mathsf{T}}
```

where the elements are the means of the electric field in each direction calculated as finite differences:

```math
\begin{aligned}
	\mathcal{E}_r^{i,j,k} &= \dfrac{1}{2}\left(\dfrac{\Phi_{i+1,j,k}-\Phi_{i,j,k}}{r_{i+1} - r_{i}} + \dfrac{\Phi_{i,j,k}-\Phi_{i-1,j,k}}{r_{i} - r_{i-1}}\right)\hspace{10pt},\\
	\mathcal{E}_{\varphi}^{i,j,k} &= \dfrac{1}{2 r_{i}}\left(\dfrac{\Phi_{i,j+1,k}-\Phi_{i,j,k}}{\varphi_{j+1} - \varphi_{j}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j-1,k}}{\varphi_{j} - \varphi_{j-1}}\right)\hspace{10pt},\\
	\mathcal{E}_z^{i,j,k} &= \dfrac{1}{2}\left( \dfrac{\Phi_{i,j,k+1}-\Phi_{i,j,k}}{z_{k+1} - z_{k}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j,k-1}}{z_{k} - z_{k-1}} \right)\hspace{10pt}.
\end{aligned}
```