	.file	"testsimd.cpp"
	.text
	.section	.rodata
.LC0:
	.string	"%f\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB4105:
	.cfi_startproc
	leaq	8(%rsp), %r10
	.cfi_def_cfa 10, 0
	andq	$-32, %rsp
	pushq	-8(%r10)
	pushq	%rbp
	.cfi_escape 0x10,0x6,0x2,0x76,0
	movq	%rsp, %rbp
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%r10
	.cfi_escape 0xf,0x3,0x76,0x58,0x6
	.cfi_escape 0x10,0xf,0x2,0x76,0x78
	.cfi_escape 0x10,0xe,0x2,0x76,0x70
	.cfi_escape 0x10,0xd,0x2,0x76,0x68
	.cfi_escape 0x10,0xc,0x2,0x76,0x60
	pushq	%rbx
	subq	$352, %rsp
	.cfi_escape 0x10,0x3,0x2,0x76,0x50
	movq	%rsp, %rax
	movq	%rax, %rbx
	movq	$1000000000, -72(%rbp)
	movq	-72(%rbp), %rax
	leaq	-1(%rax), %rdx
	movq	%rdx, -64(%rbp)
	movq	%rdx, %rax
	addq	$1, %rax
	movq	%rax, -400(%rbp)
	movq	$0, -392(%rbp)
	movq	%rdx, %rax
	addq	$1, %rax
	movq	%rax, %r14
	movl	$0, %r15d
	movq	%rdx, %rax
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %ecx
	movl	$0, %edx
	divq	%rcx
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -56(%rbp)
	movq	-72(%rbp), %rax
	subq	$1, %rax
	movq	%rax, -88(%rbp)
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r12
	movl	$0, %r13d
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r10
	movl	$0, %r11d
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %ecx
	movl	$0, %edx
	divq	%rcx
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -96(%rbp)
	movq	-72(%rbp), %rax
	subq	$1, %rax
	movq	%rax, -104(%rbp)
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %r8
	movl	$0, %r9d
	movq	%rax, %rdx
	addq	$1, %rdx
	movq	%rdx, %rsi
	movl	$0, %edi
	addq	$1, %rax
	leaq	0(,%rax,4), %rdx
	movl	$16, %eax
	subq	$1, %rax
	addq	%rdx, %rax
	movl	$16, %edi
	movl	$0, %edx
	divq	%rdi
	imulq	$16, %rax, %rax
	subq	%rax, %rsp
	movq	%rsp, %rax
	addq	$3, %rax
	shrq	$2, %rax
	salq	$2, %rax
	movq	%rax, -112(%rbp)
	movq	$0, -80(%rbp)
	movq	$8, -120(%rbp)
.L6:
	movq	-72(%rbp), %rax
	subq	-80(%rbp), %rax
	cmpq	%rax, -120(%rbp)
	ja	.L2
	movq	-80(%rbp), %rax
	leaq	0(,%rax,4), %rdx
	movq	-56(%rbp), %rax
	addq	%rdx, %rax
	movq	%rax, -384(%rbp)
	movq	-384(%rbp), %rax
	vmovups	(%rax), %xmm0
	vinsertf128	$0x1, 16(%rax), %ymm0, %ymm0
	vmovaps	%ymm0, -176(%rbp)
	movq	-80(%rbp), %rax
	leaq	0(,%rax,4), %rdx
	movq	-96(%rbp), %rax
	addq	%rdx, %rax
	movq	%rax, -376(%rbp)
	movq	-376(%rbp), %rax
	vmovups	(%rax), %xmm0
	vinsertf128	$0x1, 16(%rax), %ymm0, %ymm0
	vmovaps	%ymm0, -208(%rbp)
	vmovaps	-176(%rbp), %ymm0
	vmovaps	%ymm0, -336(%rbp)
	vmovaps	-208(%rbp), %ymm0
	vmovaps	%ymm0, -368(%rbp)
	vmovaps	-336(%rbp), %ymm0
	vaddps	-368(%rbp), %ymm0, %ymm0
	vmovaps	%ymm0, -240(%rbp)
	movq	-80(%rbp), %rax
	leaq	0(,%rax,4), %rdx
	movq	-112(%rbp), %rax
	addq	%rdx, %rax
	movq	%rax, -248(%rbp)
	vmovaps	-240(%rbp), %ymm0
	vmovaps	%ymm0, -304(%rbp)
	vmovaps	-304(%rbp), %ymm0
	movq	-248(%rbp), %rax
	vmovups	%xmm0, (%rax)
	vextractf128	$0x1, %ymm0, 16(%rax)
	nop
	movq	-120(%rbp), %rax
	addq	%rax, -80(%rbp)
	jmp	.L6
.L2:
	movq	-112(%rbp), %rax
	vmovss	(%rax), %xmm0
	vcvtss2sd	%xmm0, %xmm0, %xmm1
	vmovq	%xmm1, %rax
	vmovq	%rax, %xmm0
	leaq	.LC0(%rip), %rdi
	movl	$1, %eax
	call	printf@PLT
	movq	%rbx, %rsp
	movl	$0, %eax
	leaq	-48(%rbp), %rsp
	popq	%rbx
	popq	%r10
	.cfi_def_cfa 10, 0
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	leaq	-8(%r10), %rsp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4105:
	.size	main, .-main
	.ident	"GCC: (Debian 10.2.1-6) 10.2.1 20210110"
	.section	.note.GNU-stack,"",@progbits
