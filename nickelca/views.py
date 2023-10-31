from django.shortcuts import render

# Create your views here.

def nickelca(requests):
    return render(requests, 'nickelca.html')