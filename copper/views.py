from django.shortcuts import render

# Create your views here.

def copper(requests):
    return render(requests, 'copper.html')